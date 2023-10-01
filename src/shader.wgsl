const pi = 3.1415926535897932385;
const infinity = 1000000000.0;
const delta = 0.00000001;

const samples_per_pixel = 100;
const max_depth = 10;

const SphereNum = 4;

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
fn vec3_random() -> vec3<f32> {
    return vec3<f32>(random(), random(), random());
}
fn vec3_random_double(min: f32, max: f32) -> vec3<f32> {
    return vec3<f32>(random_double(min, max), random_double(min, max), random_double(min, max));
}

fn random_in_unit_sphere() -> vec3<f32> {
    while true {
        var v = vec3_random_double(-1.0, 1.0);
        if (length(v) < 1.0) {
            return v;
        }
    }
    return vec3<f32>();
}
fn random_unit_vector() -> vec3<f32> {
    return unit(random_in_unit_sphere());
}
fn random_on_hemisphere(normal: vec3<f32>) -> vec3<f32> {
    var on_unit_sphere = random_unit_vector();
    if (dot(on_unit_sphere, normal) > 0.0) {
        return on_unit_sphere;
    } else {
        return -on_unit_sphere;
    }
}

fn near_zero(v: vec3<f32>) -> bool {
    return length(v) < delta;
}

fn reflectance(cosine: f32, ref_idx: f32) -> f32 {
    // Use Schlick's approximation for reflectance
    var r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) * pow((1.0-cosine), 5.0);
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
    height: u32,
    width: u32,
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
fn ray_color(ray_ptr: ptr<function, Ray>, spheres: ptr<function, array<Sphere, SphereNum>>) -> vec3<f32> {
    var rec = HitRecord();
    var multiplier = vec3<f32>(1.0, 1.0, 1.0);
    for (var i = 0; i < max_depth; i = i+1) {
        if (spheres_hit(spheres, *ray_ptr, Interval(0.001, infinity), &rec)) {
            var scattered = Ray();
            var attenuation = vec3<f32>();
            if (scatter(rec.material, *ray_ptr, &rec, &attenuation, &scattered)) {
                multiplier = multiplier * attenuation;
                *ray_ptr = scattered;
                continue;
            }
            return vec3<f32>(0.0, 0.0, 0.0);
        }

        var unit_direction = unit((*ray_ptr).dir);
        var a = 0.5 * (unit_direction.z + 1.0);
        return multiplier * ((1.0 - a)*vec3<f32>(1.0, 1.0, 1.0) + a*vec3<f32>(0.5, 0.7, 1.0));
    }

    return vec3<f32>(0.0, 0.0, 0.0);
}

struct Material {
    // 0 .. lambertian, 1 .. metal, 2 .. dielectric
    material: i32,
    albedo: vec3<f32>,
    fuzz: f32,
    ir: f32,
};
fn scatter(
    material: Material,
    r_in: Ray,
    rec_ptr: ptr<function, HitRecord>,
    attenuation_ptr: ptr<function, vec3<f32>>,
    scattered_ptr: ptr<function, Ray>
) -> bool {
    switch material.material {
        case 0: { // lambertian
            var scatter_direction = (*rec_ptr).normal + random_unit_vector();

            if (near_zero(scatter_direction)) {
                scatter_direction = (*rec_ptr).normal;
            }

            *scattered_ptr = Ray((*rec_ptr).p, scatter_direction);
            *attenuation_ptr = material.albedo;
            return true;
        }
        case 1: { // metal
            var reflected = reflect(unit(r_in.dir), (*rec_ptr).normal);
            *scattered_ptr = Ray((*rec_ptr).p, reflected + material.fuzz*random_unit_vector());
            *attenuation_ptr = material.albedo;
            return dot((*scattered_ptr).dir, (*rec_ptr).normal) > 0.0;
        }
        case 2: { // dielectric
            *attenuation_ptr = vec3<f32>(1.0, 1.0, 1.0);
            var refraction_ratio = 1.0 / material.ir;
            if (!(*rec_ptr).front_face) {
                refraction_ratio = material.ir;
            }

            var unit_direction = unit(r_in.dir);
            var cos_theta = min(dot(-unit_direction, (*rec_ptr).normal), 1.0);
            var sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            var cannot_refract = refraction_ratio * sin_theta > 1.0;
            var direction = vec3<f32>();

            if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random()) {
                direction = reflect(unit_direction, (*rec_ptr).normal);
            } else {
                direction = refract(unit_direction, (*rec_ptr).normal, refraction_ratio);
            }

            *scattered_ptr = Ray((*rec_ptr).p, direction);
            return true;
        }
        default: {
            return false;
        }
    }
}

struct HitRecord {
    p: vec3<f32>,
    normal: vec3<f32>,
    material: Material,
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
    material: Material,
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
    (*rec).material = sphere.material;

    return true;
}

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

fn get_ray(position: vec2<f32>) -> Ray {
    var view_center= unit(camera.direction);
    var view_x = unit(vec3<f32>(view_center.y, -view_center.x, 0.0));
    var view_y = unit(cross(view_x, camera.direction));

    var x = position.x + random_double(-0.5, 0.5) / f32(camera.width);
    var y = position.y + random_double(-0.5, 0.5) / f32(camera.height);

    var ray_direction = view_center + view_x * x * f32(camera.width) / f32(camera.height) + view_y * y;
    return Ray(camera.position, ray_direction);
}

fn adjust_color(pixel_color: vec3<f32>, samples_per_pixel: i32) -> vec3<f32> {
    var color = pixel_color / f32(samples_per_pixel);

    return vec3<f32>(clamp(color.x, 0.0, 1.0), clamp(color.y, 0.0, 1.0), clamp(color.z, 0.0, 1.0));
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    var material_ground = Material(0, vec3<f32>(0.8, 0.8, 0.0), 0.0, 1.5);
    // var material_center = Material(0, vec3<f32>(0.7, 0.3, 0.3), 0.0, 0.0);
    var material_center = Material(2, vec3<f32>(0.7, 0.3, 0.3), 0.0, 1.5);
    // var material_left = Material(1, vec3<f32>(0.8, 0.8, 0.8), 0.3, 0.0);
    var material_left = Material(2, vec3<f32>(0.8, 0.8, 0.8), 0.0, 1.5);
    var material_right = Material(1, vec3<f32>(0.8, 0.6, 0.2), 1.0, 1.5);

    var spheres = array<Sphere, SphereNum>(
        Sphere(vec3<f32>(3.0, 0.0, 0.0), 0.5, material_center),
        Sphere(vec3<f32>(0.0, 0.0, -1000.5), 1000.0, material_ground),
        Sphere(vec3<f32>(3.0, -1.0, 0.0), 0.5, material_left),
        Sphere(vec3<f32>(3.0, 1.0, 0.0), 0.5, material_right),
    );
    seed = u32(clip(in.position.x) * 100000000.0 + clip(in.position.y) * 10000.0);

    var color = vec3<f32>(0.0, 0.0, 0.0);
    for (var sample = 0; sample < samples_per_pixel; sample = sample+1) {
        var ray = get_ray(in.position);
        color += ray_color(&ray, &spheres);
    }
    return vec4<f32>(adjust_color(color, samples_per_pixel), 1.0);
}
