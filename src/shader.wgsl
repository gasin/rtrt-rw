// Vertex shader

struct Camera {
    position: vec3<f32>,
    direction: vec3<f32>,
};
@group(0) @binding(0)
var<uniform> camera: Camera;

// const camera = Camera(vec3<f32>(0.0, 0.0, 0.0), vec3<f32>(1.0, 0.0, 0.0));

const focal_length = 1.0;
const viewpoint_height = 2.0;
// const viewpoint_width = viewpoint_height * 16.0 / 9.0;

fn hit_sphere(center: vec3<f32>, radius: f32, origin: vec3<f32>, ray: vec3<f32>) -> f32 {
    var oc = origin - center;
    var a = dot(ray, ray);
    var half_b = dot(oc, ray);
    var c = dot(oc, oc) - radius*radius;
    var discriminant = half_b*half_b - a*c;

    if (discriminant < 0.0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant) ) / a;
    }
}

struct VertexInput {
    @location(0) position: vec2<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) position: vec2<f32>,
};

@vertex
fn vs_main(
    model: VertexInput,
) -> VertexOutput {
    var out: VertexOutput;
    out.position = model.position;
    out.clip_position = vec4<f32>(model.position, 0.0, 1.0);
    return out;
}

// Fragment shader

fn clip(f: f32) -> f32 {
    return (f + 1.0) / 2.0;
}

fn clip_v3(v: vec3<f32>) -> vec3<f32> {
    return vec3<f32>(clip(v.x), clip(v.y), clip(v.z));
}

fn unit(v: vec3<f32>) -> vec3<f32> {
    return v / length(v);
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    var view_center= unit(camera.direction);
    var view_x = unit(vec3<f32>(view_center.y, -view_center.x, 0.0));
    var view_y = unit(cross(view_x, camera.direction));

    var ray_direction = view_center + view_x * in.position.x * 4.0 / 3.0 + view_y * in.position.y;

    var t = hit_sphere(vec3<f32>(3.0, 0.0, 0.0), 0.5, camera.position, ray_direction);
    if (t > 0.0) {
        return vec4<f32>(0.0, 0.0, 1.0, 1.0);
    }

    return vec4<f32>(clip_v3(unit(ray_direction)), 1.0);
}
