use std::iter;

use wgpu::util::DeviceExt;
use winit::{
    event::*,
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder}, dpi::Position,
};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Vertex {
    position: [f32; 2],
}

impl Vertex {
    fn desc() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x2,
                },
            ],
        }
    }
}

struct Camera {
    position: [f32; 3],
    direction: [f32; 3],
}

// 画面のサイズの逆比にする
const WIDTH: usize = 3;
const HEIGHT: usize = 4;

struct State {
    surface: wgpu::Surface,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    size: winit::dpi::PhysicalSize<u32>,
    render_pipeline: wgpu::RenderPipeline,
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
    window: Window,
    camera: Camera,
}

impl State {
    async fn new(window: Window) -> Self {
        let size = window.inner_size();

        // The instance is a handle to our GPU
        // BackendBit::PRIMARY => Vulkan + Metal + DX12 + Browser WebGPU
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            dx12_shader_compiler: Default::default(),
        });

        // # Safety
        //
        // The surface needs to live as long as the window that created it.
        // State owns the window so this should be safe.
        let surface = unsafe { instance.create_surface(&window) }.unwrap();

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::default(),
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .unwrap();

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: None,
                    features: wgpu::Features::empty(),
                    // WebGL doesn't support all of wgpu's features, so if
                    // we're building for the web we'll have to disable some.
                    limits: if cfg!(target_arch = "wasm32") {
                        wgpu::Limits::downlevel_webgl2_defaults()
                    } else {
                        wgpu::Limits::default()
                    },
                },
                None, // Trace path
            )
            .await
            .unwrap();

        let surface_caps = surface.get_capabilities(&adapter);
        // Shader code in this tutorial assumes an Srgb surface texture. Using a different
        // one will result all the colors comming out darker. If you want to support non
        // Srgb surfaces, you'll need to account for that when drawing to the frame.
        let surface_format = surface_caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(surface_caps.formats[0]);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width,
            height: size.height,
            present_mode: surface_caps.present_modes[0],
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
        };
        surface.configure(&device, &config);

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let render_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Render Pipeline Layout"),
                bind_group_layouts: &[],
                push_constant_ranges: &[],
            });

        let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Render Pipeline"),
            layout: Some(&render_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                buffers: &[Vertex::desc()],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState {
                        color: wgpu::BlendComponent::REPLACE,
                        alpha: wgpu::BlendComponent::REPLACE,
                    }),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                // Setting this to anything other than Fill requires Features::POLYGON_MODE_LINE
                // or Features::POLYGON_MODE_POINT
                polygon_mode: wgpu::PolygonMode::Fill,
                // Requires Features::DEPTH_CLIP_CONTROL
                unclipped_depth: false,
                // Requires Features::CONSERVATIVE_RASTERIZATION
                conservative: false,
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            // If the pipeline will be used with a multiview render pass, this
            // indicates how many array layers the attachments will have.
            multiview: None,
        });

        let mut vertices = [Vertex {
            position: [0.0, 0.0],
        }; (HEIGHT+1) * (WIDTH+1)];
        for i in 0..HEIGHT+1 {
            for j in 0..WIDTH+1 {
                let x: f32 = (i as f32) * 2.0 / (HEIGHT as f32) - 1.0;
                let y: f32 = (j as f32) * 2.0 / (WIDTH as f32) - 1.0;
                vertices[i*(WIDTH+1)+j] = Vertex {
                    position: [x, y],
                }
            }
        }

        // let mut indices: &[u16] = &[0, WIDTH as u16, 1, 1, WIDTH as u16 + 1, 2];
        let mut indices = [0; 6 * WIDTH * HEIGHT];
        for i in 0..HEIGHT {
            for j in 0..WIDTH {
                indices[6 * (i*WIDTH + j)] = (i as u16) * (WIDTH as u16 + 1) + j as u16;
                indices[6 * (i*WIDTH + j) + 1] = (i as u16 +1) * (WIDTH as u16 + 1) + j as u16;
                indices[6 * (i*WIDTH + j) + 2] = (i as u16) * (WIDTH as u16 + 1) + j as u16 + 1;
                indices[6 * (i*WIDTH + j) + 3] = (i as u16 +1) * (WIDTH as u16 + 1) + j as u16;
                indices[6 * (i*WIDTH + j) + 4] = (i as u16 +1) * (WIDTH as u16 + 1) + j as u16 + 1;
                indices[6 * (i*WIDTH + j) + 5] = (i as u16) * (WIDTH as u16 + 1) + j as u16 + 1;
            }
        }

        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Vertex Buffer"),
            // contents: bytemuck::cast_slice(VERTICES),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });
        let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Index Buffer"),
            // contents: bytemuck::cast_slice(INDICES),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        let num_indices = 6 * WIDTH as u32 * HEIGHT as u32;

        let camera = Camera{
            position: [0.0, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0],
        };

        Self {
            surface,
            device,
            queue,
            config,
            size,
            render_pipeline,
            vertex_buffer,
            index_buffer,
            num_indices,
            window,
            camera,
        }
    }

    pub fn window(&self) -> &Window {
        &self.window
    }

    pub fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        println!("{} {}", new_size.width, new_size.height);
        if new_size.width > 0 && new_size.height > 0 {
            self.size = new_size;
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
        }
    }

    #[allow(unused_variables)]
    fn input(&mut self, event: &WindowEvent) -> bool {
        false
    }

    fn update(&mut self) {}

    // fn step(&mut self) {
    //     self.step = self.step + 10;
    //     let mut vertices = [Vertex {
    //         position: [0.0, 0.0, 0.0],
    //         color: [0.0, 0.0, 0.0],
    //     }; (HEIGHT+1) * (WIDTH+1)];
    //     for i in 0..HEIGHT+1 {
    //         for j in 0..WIDTH+1 {
    //             let x: f32 = (i as f32) * 2.0 / (HEIGHT as f32) - 1.0;
    //             let y: f32 = (j as f32) * 2.0 / (WIDTH as f32) - 1.0;
    //             let pre_col_r = ((i as u32 + self.step) % (HEIGHT as u32 + 1)) as f32;
    //             let col_r = if pre_col_r*2.0 < HEIGHT as f32 { pre_col_r / HEIGHT as f32 } else { 1.0 - pre_col_r / HEIGHT as f32};
    //             let pre_col_b = ((j as u32 + self.step) % (WIDTH as u32 + 1)) as f32;
    //             let col_b = if pre_col_b*2.0 < WIDTH as f32 { pre_col_b / WIDTH as f32 } else { 1.0 - pre_col_b / WIDTH as f32};
    //             vertices[i*(WIDTH+1)+j] = Vertex {
    //                 position: [x, y, 0.0],
    //                 color: [col_r, 0.2, col_b],
    //                 // color: [i as f32 / HEIGHT as f32, 0.0, j as f32 / WIDTH as f32],
    //             };
    //             if i == 0 {
    //                 println!("{}", col_b);
    //             }
    //         }
    //     }
    //     println!("");
    //     self.queue.write_buffer(
    //         &self.vertex_buffer,
    //         0,
    //         bytemuck::cast_slice(&vertices),
    //     );
    // }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.1,
                            g: 0.2,
                            b: 0.3,
                            a: 1.0,
                        }),
                        store: true,
                    },
                })],
                depth_stencil_attachment: None,
            });

            render_pass.set_pipeline(&self.render_pipeline);
            render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
            render_pass.draw_indexed(0..self.num_indices, 0, 0..1);
        }

        self.queue.submit(iter::once(encoder.finish()));
        output.present();

        Ok(())
    }
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen(start))]
pub async fn run() {
    cfg_if::cfg_if! {
        if #[cfg(target_arch = "wasm32")] {
            std::panic::set_hook(Box::new(console_error_panic_hook::hook));
            console_log::init_with_level(log::Level::Warn).expect("Could't initialize logger");
        } else {
            env_logger::init();
        }
    }

    let event_loop = EventLoop::new();
    let window = WindowBuilder::new().build(&event_loop).unwrap();

    #[cfg(target_arch = "wasm32")]
    {
        // Winit prevents sizing with CSS, so we have to set
        // the size manually when on web.
        use winit::dpi::PhysicalSize;
        window.set_inner_size(PhysicalSize::new(450, 400));

        use winit::platform::web::WindowExtWebSys;
        web_sys::window()
            .and_then(|win| win.document())
            .and_then(|doc| {
                let dst = doc.get_element_by_id("wasm-example")?;
                let canvas = web_sys::Element::from(window.canvas());
                dst.append_child(&canvas).ok()?;
                Some(())
            })
            .expect("Couldn't append canvas to document body.");
    }

    // State::new uses async code, so we're going to wait for it to finish
    let mut state = State::new(window).await;
    // let mut last_frame_time = std::time::Instant::now();
    event_loop.run(move |event, _, control_flow| {
        match event {
            Event::WindowEvent {
                ref event,
                window_id,
            } if window_id == state.window().id() => {
                if !state.input(event) {
                    match event {
                        WindowEvent::CloseRequested
                        | WindowEvent::KeyboardInput {
                            input:
                                KeyboardInput {
                                    state: ElementState::Pressed,
                                    virtual_keycode: Some(VirtualKeyCode::Escape),
                                    ..
                                },
                            ..
                        } => *control_flow = ControlFlow::Exit,
                        WindowEvent::Resized(physical_size) => {
                            state.resize(*physical_size);
                        }
                        WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                            // new_inner_size is &mut so w have to dereference it twice
                            state.resize(**new_inner_size);
                        }
                        _ => {}
                    }
                }
            }
            Event::RedrawRequested(window_id) if window_id == state.window().id() => {
                state.update();
                match state.render() {
                    Ok(_) => {}
                    // Reconfigure the surface if it's lost or outdated
                    Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                        state.resize(state.size)
                    }
                    // The system is out of memory, we should probably quit
                    Err(wgpu::SurfaceError::OutOfMemory) => *control_flow = ControlFlow::Exit,
                    // We're ignoring timeouts
                    Err(wgpu::SurfaceError::Timeout) => log::warn!("Surface timeout"),
                }
            }
            Event::MainEventsCleared => {
                // RedrawRequested will only trigger once, unless we manually
                // request it.
                state.window().request_redraw();
            }
            _ => {
                // let current_time = std::time::Instant::now();
                // if current_time - last_frame_time >= std::time::Duration::from_millis(500) {
                //     state.step();
                //     match state.render() {
                //         Ok(_) => {}
                //         // Reconfigure the surface if it's lost or outdated
                //         Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                //             state.resize(state.size)
                //         }
                //         // The system is out of memory, we should probably quit
                //         Err(wgpu::SurfaceError::OutOfMemory) => *control_flow = ControlFlow::Exit,
                //         // We're ignoring timeouts
                //         Err(wgpu::SurfaceError::Timeout) => log::warn!("Surface timeout"),
                //     }
                //     last_frame_time = current_time;
                // }
            }
        }
    });
}