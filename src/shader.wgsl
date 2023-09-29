const pi = 3.1415926535897932385;
const infinity = 1000000000.0;

var<private> seed: u32 = 2463534242u;
fn rand_gen() -> u32 {
  seed = seed ^ (seed << 13u);
  seed = seed ^ (seed >> 17u);
  seed = seed ^ (seed << 5u);
  return seed;
}
fn random() -> f32 {
    return f32(rand_gen()) / pow(2.0, 32.0);
}
fn random_double(min: f32, max: f32) -> f32 {
    return min + (max-min) * random();
}

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

struct Interval {
    min: f32,
    max: f32,
};
fn interval_contains(interval: Interval, x: f32) -> bool {
    return interval.min <= x &&  x <= interval.max;
}
fn interval_surrounds(interval: Interval, x: f32) -> bool {
    return interval.min < x && x < interval.max;
}

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
fn ray_color(ray: Ray, spheres: ptr<function, array<Sphere, SphereNum>>) -> vec3<f32> {
    var rec = HitRecord();
    if (spheres_hit(spheres, ray, Interval(0.0, infinity), &rec)) {
        return 0.5 * (rec.normal + vec3<f32>(1.0, 1.0, 1.0));
    }

    var unit_direction = unit(ray.dir);
    var a = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - a)*vec3<f32>(1.0, 1.0, 1.0) + a*vec3<f32>(0.5, 0.7, 1.0);
}

struct HitRecord {
    p: vec3<f32>,
    normal: vec3<f32>,
    t: f32,
    front_face: bool,
};
fn hit_record_set_face_normal(hit_record_ptr: ptr<function, HitRecord>, r: Ray, outward_normal: vec3<f32>) {
    // outward_normal is assumed to have unit length

    (*hit_record_ptr).front_face = dot(r.dir, outward_normal) < 0.0;
    if ((*hit_record_ptr).front_face) {
        (*hit_record_ptr).normal = outward_normal;
    } else {
        (*hit_record_ptr).normal = -outward_normal;
    }
}

struct Sphere {
    center: vec3<f32>,
    radius: f32,
};
fn sphere_hit(sphere: Sphere, ray: Ray, ray_t: Interval, rec: ptr<function, HitRecord>) -> bool {
    var oc = ray.orig - sphere.center;
    var a = dot(ray.dir, ray.dir);
    var half_b = dot(oc, ray.dir);
    var c = dot(oc, oc) - sphere.radius*sphere.radius;
    var discriminant = half_b*half_b - a*c;
    if (discriminant < 0.0) {
        return false;
    }
    var sqrtd = sqrt(discriminant);

    var root = (-half_b - sqrtd) / a;
    if (!interval_surrounds(ray_t, root)) {
        root = (-half_b + sqrtd) / a;
        if (!interval_surrounds(ray_t, root)) {
            return false;
        }
    }

    (*rec).t = root;
    (*rec).p = ray_at(ray, (*rec).t);
    var outward_normal = ((*rec).p - sphere.center) / sphere.radius;
    hit_record_set_face_normal(rec, ray, outward_normal);

    return true;
}

const SphereNum = 2;
fn spheres_hit(spheres: ptr<function, array<Sphere, SphereNum>>, ray: Ray, ray_t: Interval, hit_record_ptr: ptr<function, HitRecord>) -> bool {
    var temp_rec = HitRecord();
    var hit_anything = false;
    var closest_so_far = ray_t.max;

    for (var i = 0; i < SphereNum; i = i+1) {
        var sphere = (*spheres)[i];
        if (sphere_hit(sphere, ray, Interval(ray_t.min, closest_so_far), &temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            *hit_record_ptr = temp_rec;
        }
    }

    return hit_anything;
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

    var spheres = array<Sphere, SphereNum>(
        Sphere(vec3<f32>(3.0, 0.0, 0.0), 0.5),
        Sphere(vec3<f32>(0.0, 0.0, -101.0), 100.0)
    );

    var color = ray_color(ray, &spheres);
    return vec4<f32>(color, 1.0);
}
