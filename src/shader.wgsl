// Vertex shader

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

struct Camera {
    position: vec3<f32>,
    direction: vec3<f32>,
};
@group(0) @binding(0)
var<uniform> camera: Camera;

struct Ray {
    orig: vec3<f32>,
    dir: vec3<f32>,
};
fn ray_at(ray: Ray, t: f32) -> vec3<f32> {
    return ray.orig + t*ray.dir;
}

struct Sphere {
    center: vec3<f32>,
    radius: f32,
};
fn sphere_hit(sphere: Sphere, ray: Ray) -> f32 {
    var oc = ray.orig - sphere.center;
    var a = dot(ray.dir, ray.dir);
    var half_b = dot(oc, ray.dir);
    var c = dot(oc, oc) - sphere.radius*sphere.radius;
    var discriminant = half_b*half_b - a*c;

    if (discriminant < 0.0) {
        return -1.0;
    } else {
        return (-half_b - sqrt(discriminant) ) / a;
    }
}

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
    var ray = Ray(camera.position, ray_direction);

    var sphere = Sphere(vec3<f32>(3.0, 0.0, 0.0), 0.5);
    var t = sphere_hit(sphere, ray);
    if (t > 0.0) {
        var n = unit(ray_at(ray, t) - sphere.center);
        return vec4<f32>(clip_v3(n), 1.0);
    }

    return vec4<f32>(clip_v3(unit(ray_direction)), 1.0);
}
