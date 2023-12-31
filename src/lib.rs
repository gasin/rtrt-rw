use std::iter;
use rand::prelude::*;

use wgpu::{util::DeviceExt, BindGroupLayoutDescriptor};
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

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Material {
    // 0 .. lambertian, 1 .. metal, 2 .. dielectric
    material: i32,
    fuzz: f32,
    _padding1: [u32; 2],
    albedo: [f32; 3],
    ir: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Sphere {
    center: [f32; 3],
    radius: f32,
    material: Material,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Camera {
    position: [f32; 3],
    _padding: u32,
    direction: [f32; 3],
    height: u32,
    width: u32,
    _padding2: u32,
    _padding3: u32,
    _padding4: u32,
}

struct CameraController {
    speed: f32,
    is_forward_pressed: bool,
    is_backward_pressed: bool,
    is_left_pressed: bool,
    is_right_pressed: bool,
    is_up_pressed: bool,
    is_down_pressed: bool,
}

impl CameraController {
    fn new(speed: f32) -> Self {
        Self {
            speed,
            is_forward_pressed: false,
            is_backward_pressed: false,
            is_left_pressed: false,
            is_right_pressed: false,
            is_up_pressed: false,
            is_down_pressed: false,
        }
    }

    fn process_events(&mut self, event: &WindowEvent) -> bool {
        match event {
            WindowEvent::KeyboardInput {
                input: KeyboardInput {
                    state,
                    virtual_keycode: Some(keycode),
                    ..
                },
                ..
            } => {
                let is_pressed = *state == ElementState::Pressed;
                match keycode {
                    VirtualKeyCode::W | VirtualKeyCode::Up => {
                        self.is_forward_pressed = is_pressed;
                        true
                    }
                    VirtualKeyCode::A | VirtualKeyCode::Left => {
                        self.is_left_pressed = is_pressed;
                        true
                    }
                    VirtualKeyCode::S | VirtualKeyCode::Down => {
                        self.is_backward_pressed = is_pressed;
                        true
                    }
                    VirtualKeyCode::D | VirtualKeyCode::Right => {
                        self.is_right_pressed = is_pressed;
                        true
                    }
                    VirtualKeyCode::F => {
                        self.is_up_pressed = is_pressed;
                        true
                    }
                    VirtualKeyCode::G => {
                        self.is_down_pressed = is_pressed;
                        true
                    }
                    _ => false,
                }
            }
            _ => false,
        }
    }

    fn update_camera(&self, camera: &mut Camera) {

        let direction = [
            camera.position[0] - camera.direction[0],
            camera.position[1] - camera.direction[1],
            camera.position[2] - camera.direction[2],
        ];

        if self.is_forward_pressed {
            camera.position[0] += direction[0] * self.speed;
            camera.position[1] += direction[1] * self.speed;
        }
        if self.is_backward_pressed {
            camera.position[0] -= direction[0] * self.speed;
            camera.position[1] -= direction[1] * self.speed;
        }
        if self.is_right_pressed {
            camera.position[0] += direction[1] * self.speed;
            camera.position[1] -= direction[0] * self.speed;
        }
        if self.is_left_pressed {
            camera.position[0] -= direction[1] * self.speed;
            camera.position[1] += direction[0] * self.speed;
        }
        if self.is_up_pressed {
            camera.position[2] += self.speed;
        }
        if self.is_down_pressed {
            camera.position[2] -= self.speed;
        }
    }
}

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
    sphere_bind_group: wgpu::BindGroup,
    camera: Camera,
    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,
    camera_controller: CameraController,
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

        let vertices = [
            Vertex { position: [-1.0,-1.0] },
            Vertex { position: [-1.0, 1.0] },
            Vertex { position: [ 1.0,-1.0] },
            Vertex { position: [ 1.0, 1.0] },
        ];
        let indices: [u16; 6] = [0, 2, 1, 2, 3, 1];
        let num_indices: u32 = 6;

        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Vertex Buffer"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });
        let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Index Buffer"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        let mut spheres: Vec<Sphere> = Vec::new();
        let mut rng = rand::thread_rng();
        for i in -2..=2 {
            for j in -2..=2 {
                let choose_mat: f32 = rng.gen();
                let center = [(i as f32) *2.0 + 0.9*rng.gen::<f32>(), (j as f32)*2.0 + 0.9*rng.gen::<f32>(), 0.2];

                if choose_mat < 0.8 {
                    // diffuse
                    let albedo = [
                        rng.gen::<f32>()*rng.gen::<f32>(),
                        rng.gen::<f32>()*rng.gen::<f32>(),
                        rng.gen::<f32>()*rng.gen::<f32>(),
                    ];
                    let sphere_material = Material {
                        material: 0,
                        fuzz: 0.0,
                        _padding1: [0, 0],
                        albedo: albedo,
                        ir: 0.0,
                    };
                    spheres.push(Sphere {
                        center: center,
                        radius: 0.2,
                        material: sphere_material,
                    });
                } else if choose_mat < 0.95 {
                    // metal
                    let albedo = [rng.gen::<f32>() / 2.0 + 0.5, rng.gen::<f32>() / 2.0 + 0.5, rng.gen::<f32>() / 2.0 + 0.5];
                    let fuzz = rng.gen::<f32>() / 2.0;
                    let sphere_material = Material {
                        material: 1,
                        fuzz: fuzz,
                        _padding1: [0, 0],
                        albedo: albedo,
                        ir: 0.0,
                    };
                    spheres.push(Sphere {
                        center: center,
                        radius: 0.2,
                        material: sphere_material,
                    });
                } else {
                    // glass
                    let sphere_material = Material {
                        material: 2,
                        fuzz: 0.0,
                        _padding1: [0, 0],
                        albedo: [0.0, 0.0, 0.0],
                        ir: 1.5,
                    };
                    spheres.push(Sphere {
                        center: center,
                        radius: 0.2,
                        material: sphere_material,
                    });
                }
            }
        }
        spheres.push(
            Sphere {
                center: [0.0, 0.0, -1000.0],
                radius: 1000.0,
                material: Material{
                    material: 0,
                    fuzz: 0.0,
                    _padding1: [0, 0],
                    albedo: [0.5, 0.5, 0.5],
                    ir: 1.5,
                },
            },
        );
        spheres.push(
            Sphere {
                center: [0.0, 0.0, 1.0],
                radius: 1.0,
                material: Material{
                    material: 2,
                    fuzz: 0.0,
                    _padding1: [0, 0],
                    albedo: [0.5, 0.5, 0.5],
                    ir: 1.5,
                }
            },
        );
        spheres.push(
            Sphere {
                center: [4.0, 0.0, 1.0],
                radius: 1.0,
                material: Material{
                    material: 0,
                    fuzz: 0.0,
                    _padding1: [0, 0],
                    albedo: [0.4, 0.2, 0.1],
                    ir: 0.0,
                },
            },
        );
        spheres.push(
            Sphere {
                center: [-4.0, 0.0, 1.0],
                radius: 1.0,
                material: Material{
                    material: 1,
                    fuzz: 0.0,
                    _padding1: [0, 0],
                    albedo: [0.7, 0.6, 0.5],
                    ir: 0.0,
                },
            }
        );
        let sphere_buffer = device.create_buffer_init(
            &wgpu::util::BufferInitDescriptor {
                label: Some("Material Buffer"),
                contents: bytemuck::cast_slice(&spheres),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            }
        );
        let sphere_bind_group_layout = device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage {
                            read_only: true,
                        },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }
            ],
            label: Some("sphere_bind_group_layout"),
        });
        let sphere_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            layout: &sphere_bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: sphere_buffer.as_entire_binding(),
                }
            ],
            label: Some("sphere_bind_group"),
        });

        let camera = Camera{
            position: [6.0, 5.0, 2.0],
            _padding: 0,
            direction: [0.0, 0.0, 0.0],
            height: size.height,
            width: size.width,
            _padding2: 0,
            _padding3: 0,
            _padding4: 0,
        };
        let camera_buffer = device.create_buffer_init(
            &wgpu::util::BufferInitDescriptor {
                label: Some("Camera Buffer"),
                contents: bytemuck::cast_slice(&[camera]),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            }
        );
        let camera_bind_group_layout = device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }
            ],
            label: Some("camera_bind_group_layout"),
        });
        let camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            layout: &camera_bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: camera_buffer.as_entire_binding(),
                }
            ],
            label: Some("camera_bind_group"),
        });

        let camera_controller = CameraController::new(0.2);

        let render_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Render Pipeline Layout"),
                bind_group_layouts: &[&camera_bind_group_layout, &sphere_bind_group_layout],
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
            sphere_bind_group,
            camera,
            camera_buffer,
            camera_bind_group,
            camera_controller,
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

            self.camera.height = new_size.height;
            self.camera.width = new_size.width;
            self.queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[self.camera]));
        }
    }

    #[allow(unused_variables)]
    fn input(&mut self, event: &WindowEvent) -> bool {
        self.camera_controller.process_events(event)
    }

    fn update(&mut self) {
        self.camera_controller.update_camera(&mut self.camera);
        self.queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[self.camera]));
    }

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
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.sphere_bind_group, &[]);
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
