#include <GLUT/glut.h>
#include <OPENGL/gl.h>
#include <stdio.h>
#include <math.h>
#include <vector>

//smallest cof : 이걸로 # of bounces in tracing a ray 조절가능 :D
float MIN = 0.01;

using namespace std;

float point_size;

//Resolution
float width_inc, height_inc, theta_inc, PI_inc;

//View volume
float eyex, eyey, eyez, lookAt_x, lookAt_y, lookAt_z, up_x, up_y, up_z;
float near_plane, far_plane;
float width, height;

//shading
float ambient, diffuse, specular;

//reflect & refract RGB
float refl_R, refl_G, refl_B;
float refr_R, refr_G, refr_B;

//lineage write RGB
float rec_refl_R, rec_refl_G, rec_refl_B;
float rec_refr_R, rec_refr_G, rec_refr_B;

bool is_rec_calculating;
bool is_reflected, is_refracted, is_rec_reflected, is_rec_refracted, in_shadow;

//point light attenuation
float constant_atten, linear_atten, quadratic_atten;

//depth
bool behind = false, right = true;

enum COLOR
{
    RED,
    GREEN,
    BLUE
};

enum SHAPE_TYPE
{
    SPHERE,
    TRIANGLE
};

enum INTERSECTION_TYPE
{
    REFLECTION,
    REFRACTION
};


//ray tracing을 위해 어떤 물체와 충돌했는지 번호 기록
//sphere인지 triangle인지도 기록 -> 굴절,반사체크가 다르니까
int which_recursive_number;
SHAPE_TYPE which_recursive_shape;


//glm.h 에서 가져온 코드 -> 계산 편하게 하려고 ㅎㅎ
// vector3 or position vector
class vector3
{
public:
    float x,y,z;
    
    vector3(){}
    
    vector3(float x, float y, float z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    vector3 operator-( vector3 v )
    {
        vector3 w;
        w.x = x - v.x;
        w.y = y - v.y;
        w.z = z - v.z;
        return w;
    }
    
    vector3 operator+( vector3 v )
    {
        vector3 w;
        w.x = x + v.x;
        w.y = y + v.y;
        w.z = z + v.z;
        return w;
    }
    
    // scalar product or dot product
    float operator*( vector3 &v )
    {
        return x*v.x + y*v.y + z*v.z;
    }
    
    // multiplication with scalar
    vector3 operator *( float i )
    {
        vector3 temp;
        temp.x = x*i;
        temp.y = y*i;
        temp.z = z*i;
        return temp;
    }
    
};

//direction은 vec3_unit로
class light
{
public:
    vector3 position;
    float intensity;
    light()
    {
        
    }
    light (vector3 position, float intensity)
    {
        this->position = position;
        this->intensity = intensity;
    }
};

//distance between 2 vertices
float distance(vector3 v1, vector3 v2)
{
    return sqrt( pow((v1.x-v2.x),2) +pow((v1.y-v2.y),2) +pow((v1.z-v2.z),2) );
}

//시작점 끝점을 넣으면 그 방향벡터 && 단위벡터를 꺼내줌
class vec3_unit
{
public:
    float x, y, z;
    
    vec3_unit(){}
    
    //line from u to v is the direction and we are computing unit vector in the same direction
    vec3_unit(vector3 u, vector3 v)
    {
        float dist = distance(u,v);
        
        x = v.x - u.x;
        y = v.y - u.y;
        z = v.z - u.z;
        
        x = x / dist;
        y = y / dist;
        z = z / dist;
    }
    
    // direction 구하는 vector -> 바로 단위벡터로
    vec3_unit( vector3 *u, vector3 *v)
    {
        float dist = distance( *u, *v);
        
        x = v->x - u->x;
        y = v->y - u->y;
        z = v->z - u->z;
        
        x = x / dist;
        y = y / dist;
        z = z / dist;
    }
    
    // u -> v is the direction, we are dividing it with its length to get unit vector
    // v의 unit vector 구하기
    void unit_vectorize()
    {
        vector3 u = vector3(0.0 ,0.0 ,0.0);
        vector3 v = vector3(x, y, z);
        float dist = distance( u, v);
        x /= dist;
        y /= dist;
        z /= dist;
    }
    
    //from glm.h...클래스 다르면 operator 같이 못쓰는거 귀찮! ㅜㅜ
    vec3_unit operator-( vec3_unit &v )
    {
        vec3_unit w;
        w.x = x - v.x;
        w.y = y - v.y;
        w.z = z - v.z;
        return w;
    }
    
    vec3_unit operator+( vec3_unit &v )
    {
        vec3_unit w;
        w.x = x + v.x;
        w.y = y + v.y;
        w.z = z + v.z;
        return w;
    }
    
    //이 operation에서 바로 cos값 나옴
    vec3_unit operator*( float i )
    {
        vec3_unit w;
        w.x = x * i;
        w.y = y * i;
        w.z = z * i;
        return w;
    }
    
    // inverse the direction
    void operator-()
    {
        x = -x;
        y = -y;
        z = -z;
    }
    
    float check_mag()
    {
        return x*x + y*y + z*z;
    }
};

bool check_intersection ( vector3 surface_position, SHAPE_TYPE s, int i, vec3_unit transmit_dir, enum INTERSECTION_TYPE in);

class sphere
{
public:
    float radius;
    vector3 center;
    
    //color rgb
    float cr, cg, cb;
    //shade K
    float ka, kd, ks, n;
    //반사 상수(0이면 안하는거) & index
    float refl_cof, refr_cof, refr_ind;
    
    sphere()
    {
        
    }
    //초기화
    sphere(float radius, vector3 center, float cr, float cg, float cb, float ka, float kd, float ks, float n, float refl_cof, float refr_cof, float refr_ind )
    {
        this->radius = radius;
        this->center = center;
        this->cr = cr;
        this->cg = cg;
        this->cb = cb;
        this->ka = ka;
        this->kd = kd;
        this->ks = ks;
        this->n = n;
        this->refl_cof = refl_cof;
        this->refr_cof = refr_cof;
        this->refr_ind = refr_ind;
    }
};


class triangle
{
public:
    vector3 v1, v2, v3;
    //normal vector로 구해지는 평면의 방정식 계수
    float eqa, eqb, eqc, eqd;
    vec3_unit surface_normal;
    
    float cr, cg, cb;
    float ka, kd, ks, n;
    float refl_cof, refr_cof, refr_ind;
    
    triangle(){}
    
    triangle(float eqa, float eqb, float eqc, float eqd, vec3_unit surface_normal, vector3 *v1, vector3 *v2, vector3 *v3, float cr, float cg, float cb, float ka, float kd, float ks, float n, float refl_cof, float refr_cof, float refr_ind)
    {
        this->eqa = eqa;
        this->eqb = eqb;
        this->eqc = eqc;
        this->eqd = eqd;
        this->surface_normal = surface_normal;
        this->v1 = *v1;
        this->v2 = *v2;
        this->v3 = *v3;
        this->cr = cr;
        this->cg = cg;
        this->cb = cb;
        this->ka = ka;
        this->kd = kd;
        this->ks = ks;
        this->n = n;
        this->refl_cof = refl_cof;
        this->refr_cof = refr_cof;
        this->refr_ind = refr_ind;
    }
};

// u.v cos (theta) where theta is the angle between 2 directions
float dot_product( vec3_unit u, vec3_unit v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}
//normal vector 계산할때 쓰려고
vec3_unit cross_product( vec3_unit u, vec3_unit v )
{
    vec3_unit w;
    w.x = u.y*v.z - u.z*v.y;
    w.y = u.z*v.x - u.x*v.z;
    w.z = u.x*v.y - u.y*v.x;
    return w;
}

vector3 cross_product_vector3(vector3 u, vector3 v)
{
    vector3 w;
    w.x = u.y*v.z - u.z*v.y;
    w.y = u.z*v.x - u.x*v.z;
    w.z = u.x*v.y - u.y*v.x;
    return w;
    
}

//edit on file
int no_of_spheres, no_of_triangles, no_of_light_source;
//읽어서 여기다 저장.....선언을 하면 뭐하나 쓰지를 못하는데ㅜㅜ
vector<sphere> spheres;
vector<triangle> triangles;
vector<light> lights;

//안보이면 return false
bool is_visible(vector3 surface_position, SHAPE_TYPE s, int i)
{
    vector3 *eye = new vector3 (eyex, eyey, eyez);
    
    //surface position : 현재 충돌 확인 하려고 하는 vector3 -> 즉 이것이 ray의 방향 벡터
    vec3_unit *dir = new vec3_unit( surface_position, *eye);
    
    vector3 start = surface_position;
    vector3 *direction = new vector3 (dir->x, dir->y, dir->z);
    
    float t, eye_t;
    
    //this is the 'real' direction of ray from surface to eye
    
    //P = P0 + V*t
    vector3 ray = start + *direction * t;
    
    //any object with t < eye_t is an obstacle and will make this surface under non visible
    //맨 처음 t값. t = (ray - start)/direction인데 맨 처음 ray는 eye 방향으로 나오기 떄문에
    eye_t = (*eye - start).x / direction->x;
    
    //checking for each sphere
    for( int j = 0; j < no_of_spheres; j++)
    {
        if ( s == SPHERE && i == j) continue;
        
        float local_t;
        
        /*
         v = dot_product( EO, V )
         disc = r2 - ((dot_product( EO, EO ) - v2 )
         if ( disc < 0 )
         no intersection
         else
         d = sqrt( disc )
         P = E + (v - d) * V
         */
        
        //충돌체크크크
        float A = pow(dir->x, 2) + pow (dir->y, 2) + pow (dir->z, 2);
        float B = 2 * (dir->x*(start.x - spheres[j].center.x) + dir->y*(start.y - spheres[j].center.y) + dir->z*(start.z - spheres[j].center.z));
        float C = pow(start.x - spheres[j].center.x, 2) + pow(start.y - spheres[j].center.y, 2) + pow(start.z - spheres[j].center.z, 2);
        
        //판별식
        float D = B*B - 4*A*C;
        
        if (D < 0) continue;
        float t0 = (-B - sqrt(D)) / 2*A;
        float t1 = (-B + sqrt(D)) / 2*A;
        if (t0 < t1) local_t = t0; else local_t = t1;
        
        //물체가 eye와 surface 사이에 존재 -> 안보여!
        if( local_t > 0 && local_t < eye_t) return false;
        
    }
    
    //checking for each triangle
    for (int j=0; j < no_of_triangles; j++)
    {
        if(s == TRIANGLE && i == j) continue;
        
        float local_t;
        
        //normal vector로 구해지는 평면의 방정식 계수들
        vec3_unit pn = triangles[j].surface_normal;
        float eqa = triangles[j].eqa;
        float eqb = triangles[j].eqb;
        float eqc = triangles[j].eqc;
        float eqd = triangles[j].eqd;
        
        //평면과 직선이 일치하면 pass
        if( eqa*dir->x + eqb*dir->y + eqc*dir->z == 0) continue;
        
        //평면의 방정식에 start값 대입
        local_t = eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        //기본값 대입 -> t 구하기 위해
        local_t /= (eqa*dir->x + eqb*dir->y + eqc*dir->z);
        local_t = -local_t;
        
        vec3_unit surface_normal = pn;
        
        // P = P0 + V*t
        ray = start + *direction * local_t;
        
        //containment test
        vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
        
        vector3 cp1 = cross_product_vector3(v1-v2, v3-v2);
        vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
        
        if( cp1*cp2 < 0 ) continue;
        
        vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
        vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
        
        if( cp3*cp4 < 0 ) continue;
        
        vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
        vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
        
        if( cp5*cp6 < 0 ) continue;
        
        //물체가 eye와 surface 사이에 존재 -> 안보여!
        if( local_t > 0 && local_t < eye_t ) return false;
    }
    
    return true;
    
}


//위에꼐 눈 기준으로 만나고 안만나고(보이고 안보이고)를 결정했다면, 이 함수는 빛 기준으로 닿고 안닿고를 결정
//아까는 만나면 보이는게 true, 이번엔 안만나야 그림자가 true;
bool is_shadow(vector3 surface_position, vector3 light_source_position, SHAPE_TYPE s, int i)
{
    //unit vector from surface to light source. eye에서 light로만 바뀌고 위랑 복ㅋ붙ㅋ
    vec3_unit *dir = new vec3_unit(surface_position, light_source_position);
    vector3 start = surface_position;
    vector3 *direction = new vector3 (dir->x, dir->y, dir->z);
    
    float t, light_t;
    
    //빛의 직선방정식 P = P0 + Vt
    vector3 ray = start + *direction * t;
    
    //t < light_t 인 오브젝트들이 가리니까 shadow를 만든다
    light_t = (light_source_position - start).x / direction->x;
    
    //checking for each sphere
    for( int j = 0; j < no_of_spheres; j++ )
    {
        if( s == SPHERE && i == j ) continue;
        
        float local_t;
        
        float A = pow( dir->x, 2) + pow( dir->y, 2) + pow( dir->z, 2);
        float B = 2 * ( dir->x*(start.x - spheres[j].center.x) + dir->y*(start.y - spheres[j].center.y) + dir->z*(start.z - spheres[j].center.z) );
        float C = pow( start.x - spheres[j].center.x, 2) + pow( start.y - spheres[j].center.y, 2) + pow( start.z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
        float D = B*B - 4*A*C;
        
        if( D < 0 ) continue;
        
        float t0 = (-B - sqrt(D))/2*A;
        float t1 = (-B + sqrt(D))/2*A;
        if( t0 < t1 ) local_t = t0; else local_t = t1;
        
        //물체가 eye와 surface 사이에 존재 -> 그림자
        if( local_t > 0 && local_t < light_t ) return true;
        
    }
    
    // checking for each triangle
    for( int j = 0; j < no_of_triangles; j++ )
    {
        if( s == TRIANGLE && i == j ) continue;
        
        float local_t;
        
        vec3_unit pn = triangles[j].surface_normal;
        float eqa = triangles[j].eqa;
        float eqb = triangles[j].eqb;
        float eqc = triangles[j].eqc;
        float eqd = triangles[j].eqd;
        
        if( eqa*dir->x + eqb*dir->y + eqc*dir->z == 0 ) continue;
        
        local_t =  eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        local_t /= eqa*dir->x + eqb*dir->y + eqc*dir->z;
        local_t = -local_t;
        
        vec3_unit surface_normal = pn;
        
        ray = start + *direction * local_t;
        
        // containment test
        vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
        
        vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
        vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
        
        if( cp1*cp2 < 0 ) continue;
        
        vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
        vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
        
        if( cp3*cp4 < 0 ) continue;
        
        vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
        vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
        
        if( cp5*cp6 < 0 ) continue;
        
        //물체가 eye와 surface 사이에 존재 -> 안보여!
        if( local_t > 0 && local_t < light_t ) return true;
    }
    
    return false;
}

//recursive 함수
float calculate_light_intensity( COLOR c, SHAPE_TYPE s, int i, vector3 surface_position, vec3_unit surface_normal, vec3_unit view_dir, vec3_unit transmitted_direction)
{
    //맨 처음에 만나면
    if( !is_rec_calculating )
    {
        is_rec_calculating = true;
        
        //벡터 계산 하려고 뒤집뒤집
        vec3_unit inverse_transmitted_direction = transmitted_direction;
        -inverse_transmitted_direction;
        
        if(s == TRIANGLE)
        {
            //굴절광은 그대로 들어가요 -> 얇은 평면이니까
            vec3_unit refract_dir = transmitted_direction;
            if( triangles[i].refr_cof > MIN) check_intersection(surface_position, s, i, refract_dir, REFRACTION);
            
            //this is the reflection direction
            vec3_unit reflect_dir;
            
            //반사광 구하기
            vec3_unit temp = surface_normal * (2 *  dot_product(inverse_transmitted_direction, surface_normal));
            reflect_dir = temp - inverse_transmitted_direction;
            
            //재귀함수 콜>_<
            if ( triangles[i].refl_cof > MIN ) check_intersection(surface_position, s, i, reflect_dir, REFLECTION);
            
        }
        
        
        if( s == SPHERE )
        {
            
            //REFRACTION THOROUGH SPHERE
            //to find direction of refracted ray
            
            //nt * sin(t) = ni * sin(i)
            float cos_I = dot_product( inverse_transmitted_direction, surface_normal );
            float sin_I = 1 - pow(cos_I, 2);
            float refr_index = 1/(spheres[i].refr_ind);
            float sin_T = refr_index * sin_I;
            float cos_T = 1 - pow(sin_T, 2);
            
            vec3_unit temp_I = transmitted_direction;
            
            vector3 *I = new vector3 (temp_I.x, temp_I.y, temp_I.z);
            vector3 *N = new vector3 (surface_normal.x, surface_normal.y, surface_normal.z);
            
            //this is the refraction direction from surface to in the direction of refracted ray
            vector3 *direction = new vector3();
            
            //Rr = (n * I) + (n * cosi - cost) * N, N = refr_index
            *direction = (*I + *N * cos_I) * refr_index - *N * cos_T;
            vec3_unit *refract_dir = new vec3_unit( new vector3(0.0,0.0,0.0), direction);
            
            //t is at infinity -> 원에서 만나는 점을 찾기 위해
            float t = 1000;
            vector3 start = surface_position;
            direction = new vector3( refract_dir->x, refract_dir->y, refract_dir->z);
            
            //구 안으로 들어가는 ray
            vector3 ray = start + *direction * t;
            
            //충돌지점의 t
            float local_t;
            
            //충돌체크
            float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
            float B = 2 * ( refract_dir->x*(start.x - spheres[i].center.x) + refract_dir->y*(start.y - spheres[i].center.y) + refract_dir->z*(start.z - spheres[i].center.z) );
            float C = pow( start.x - spheres[i].center.x, 2) + pow( start.y - spheres[i].center.y, 2) + pow( start.z - spheres[i].center.z, 2) - pow( spheres[i].radius, 2);
            float D = B*B - 4*A*C;
            
            float t0 = (-B - sqrt(D))/2*A;
            float t1 = (-B + sqrt(D))/2*A;
            if( t0 < t1 ) local_t = t1; else local_t = t0;
            
            t = local_t;
            
            ray = start + *direction * t;
            //원의 중심
            vector3 *cen = new vector3( spheres[i].center.x, spheres[i].center.y, spheres[i].center.z);
            //원의 법선벡터
            vec3_unit *surface_normal_1 = new vec3_unit(ray, *cen);
            //빛->표면
            vec3_unit *view_dir_1 = new vec3_unit(ray, surface_position);
            
            //to fined the direction of refracted ray
            cos_I    = dot_product( *view_dir_1, *surface_normal_1);
            sin_I    = sqrt(1 - pow(cos_I, 2));
            refr_index = spheres[i].refr_ind;
            sin_T    = refr_index * sin_I;
            cos_T    = sqrt(1 - pow(sin_T, 2));
            
            //방향 맞추려고
            temp_I = *view_dir_1;
            -temp_I;
            
            /*
             n1 = index of refraction of original medium
             n2 = index of refraction of new medium
             n = n1 / n2
             c2 = sqrt( 1 - n2 * (1 - c12) )
             
             Rr = (n * V) + (n * c1 - c2) * N
             */
            
            I = new vector3(temp_I.x, temp_I.y, temp_I.z);
            N = new vector3( surface_normal_1->x, surface_normal_1->y, surface_normal_1->z);
            
            *direction = (*I + *N * cos_I) * refr_index - *N * cos_T;
            
            refract_dir = new vec3_unit( new vector3(0,0,0), direction);
            
            if(spheres[i].refr_cof > MIN ) check_intersection(surface_position, s, i, *refract_dir, REFRACTION);
            
            //reflect direction
            vec3_unit reflect_dir;
            
            /*
             c1 = -dot_product( N, V )
             Rl = V + (2 * N * c1 )
             */
            vec3_unit temp = surface_normal * (2 * dot_product( inverse_transmitted_direction, surface_normal));
            reflect_dir = temp - inverse_transmitted_direction;
            
            if(spheres[i].refl_cof > MIN ) check_intersection(surface_position, s, i, reflect_dir, REFLECTION);
            
        }
        
        is_rec_calculating = false;
    }
    
    //local light source
    //빛마다 diffuse, specular값 누적
    float temp_diffuse = 0.0, temp_specular = 0.0;
    for( int j = 0; j < no_of_light_source; j++)
    {
        //좌표의 diffuse and specular
        float local_diffuse = 0.0, local_specular = 0.0;
        
        //light source position
        vector3 *light_source_position = new vector3 (lights[j].position.x, lights[j].position.y, lights[j].position.z);
        
        //direction : 반대로(계산)
        vec3_unit *light_dir = new vec3_unit(*light_source_position, surface_position);
        
        -(*light_dir);
        
        //그림자 check
        if (is_shadow(surface_position, *light_source_position, s, i)) continue;
        
        //diffuse -> kd * Light
        if ( s== SPHERE)
        {
            float diffuse_dot_prod = dot_product (*light_dir, surface_normal);
            
            //빛과 노말벡터가 90도 이상이면 -> 안만남 -> 안보임
            if (diffuse_dot_prod < 0) continue;
            local_diffuse = spheres[i].kd * lights[j].intensity * diffuse_dot_prod;
        }
        else // triangle
        {
            float diffuse_dot_prod;
            
            //view_dir : ray -> eye 시선방향
            if(dot_product(view_dir, surface_normal) > 0 && dot_product(*light_dir, surface_normal) > 0)
            {
                diffuse_dot_prod = dot_product (*light_dir, surface_normal);
                local_diffuse = triangles[i].kd * lights[j].intensity * diffuse_dot_prod;
            }
            
            else if ( dot_product(view_dir, surface_normal) < 0 && dot_product (*light_dir, surface_normal) < 0)
            {
                diffuse_dot_prod = -dot_product(*light_dir, surface_normal);
                local_diffuse = triangles[i].kd * lights[j].intensity * diffuse_dot_prod;
            }
            else
                continue;
        }
        
        vec3_unit half_way_vector = *light_dir + view_dir;
        half_way_vector.unit_vectorize();
        
        if( s == SPHERE )
        {
            float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), spheres[i].n);
            local_specular += spheres[i].ks * lights[j].intensity * specular_dot_prod;
        }
        else
        {
            float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), triangles[i].n);
            local_specular += triangles[i].ks * lights[j].intensity * specular_dot_prod;
        }
        
        //attenuation 적용
        local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, surface_position) + (quadratic_atten)*distance(*light_source_position, surface_position)*distance(*light_source_position, surface_position);
        local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, surface_position) + (quadratic_atten)*distance(*light_source_position, surface_position)*distance(*light_source_position, surface_position);
        
        temp_diffuse += local_diffuse;
        temp_specular += local_specular;
    }
    
    if( is_rec_calculating )
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cr + temp_specular;
            if( c == GREEN ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cg + temp_specular;
            if( c == BLUE  ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cb + temp_specular;
        }
        else
        {
            if( c == RED   ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cr + temp_specular;
            if( c == GREEN ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cg + temp_specular;
            if( c == BLUE  ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cb + temp_specular;
        }
    }
    
    if( is_rec_reflected )
    {
        if( is_rec_refracted )
        {
            if( s == SPHERE )
            {
                if( c == RED   ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cr * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + temp_specular + rec_refl_R   * ( spheres[i].refl_cof ) + rec_refr_R   * ( spheres[i].refr_cof );
                if( c == GREEN ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cg * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + temp_specular + rec_refl_G * ( spheres[i].refl_cof ) + rec_refr_G * ( spheres[i].refr_cof );
                if( c == BLUE  ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cb * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + temp_specular + rec_refl_B  * ( spheres[i].refl_cof ) + rec_refr_B  * ( spheres[i].refr_cof );
            }
            else
            {
                if( c == RED   ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cr * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + temp_specular + rec_refl_R   * ( triangles[i].refl_cof ) + rec_refr_R   * ( triangles[i].refr_cof );
                if( c == GREEN ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cg * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + temp_specular + rec_refl_G * ( triangles[i].refl_cof ) + rec_refr_G * ( triangles[i].refr_cof );
                if( c == BLUE  ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cb * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + temp_specular + rec_refl_B  * ( triangles[i].refl_cof ) + rec_refr_B  * ( triangles[i].refr_cof );
                
            }
        }
        else
        {
            
            if( s == SPHERE )
            {
                if( c == RED   ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cr * (1 - spheres[i].refl_cof ) + temp_specular + rec_refl_R   * ( spheres[i].refl_cof );
                if( c == GREEN ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cg * (1 - spheres[i].refl_cof ) + temp_specular + rec_refl_G * ( spheres[i].refl_cof );
                if( c == BLUE  ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cb * (1 - spheres[i].refl_cof ) + temp_specular + rec_refl_B  * ( spheres[i].refl_cof );
            }
            else
            {
                if( c == RED   ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cr * (1 - triangles[i].refl_cof ) + temp_specular + rec_refl_R   * ( triangles[i].refl_cof );
                if( c == GREEN ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cg * (1 - triangles[i].refl_cof ) + temp_specular + rec_refl_G * ( triangles[i].refl_cof );
                if( c == BLUE  ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cb * (1 - triangles[i].refl_cof ) + temp_specular + rec_refl_B  * ( triangles[i].refl_cof );
            }
        }
    }
    else if( is_rec_refracted )
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cr * (1 - spheres[i].refr_cof ) + temp_specular + rec_refr_R   * ( spheres[i].refr_cof );
            if( c == GREEN ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cg * (1 - spheres[i].refr_cof ) + temp_specular + rec_refr_G * ( spheres[i].refr_cof );
            if( c == BLUE  ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cb * (1 - spheres[i].refr_cof ) + temp_specular + rec_refr_B  * ( spheres[i].refr_cof );
        }
        else
        {
            if( c == RED   ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cr * (1 - triangles[i].refr_cof ) + temp_specular + rec_refr_R   * ( triangles[i].refr_cof );
            if( c == GREEN ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cg * (1 - triangles[i].refr_cof ) + temp_specular + rec_refr_G * ( triangles[i].refr_cof );
            if( c == BLUE  ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cb * (1 - triangles[i].refr_cof ) + temp_specular + rec_refr_B  * ( triangles[i].refr_cof );
        }
    }
    else
    {
        if( s == SPHERE )
        {
            if( c == RED   ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cr + temp_specular;
            if( c == GREEN ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cg + temp_specular;
            if( c == BLUE  ) return (ambient * spheres[i].ka + temp_diffuse) * spheres[i].cb + temp_specular;
        }
        else
        {
            if( c == RED   ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cr + temp_specular;
            if( c == GREEN ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cg + temp_specular;
            if( c == BLUE  ) return (ambient * triangles[i].ka + temp_diffuse) * triangles[i].cb + temp_specular;
        }
    }
    
    return ambient;
}

bool check_intersection( vector3 surface_position, SHAPE_TYPE s, int i, vec3_unit transmit_dir, enum INTERSECTION_TYPE in )
{
    
    vector3 start = surface_position;
    vector3 *eye_position = new vector3(eyex, eyey, eyez);
    
    vector3 *direction = new vector3( transmit_dir.x, transmit_dir.y, transmit_dir.z );
    
    float t = 1000;
    vector3 ray = start + *direction * t;
    
    for( int j = 0; j < no_of_spheres; j++ )
    {
        if( s == SPHERE && i == j ) continue;
        
        float local_t;
        
        float A = pow( transmit_dir.x, 2) + pow( transmit_dir.y, 2) + pow( transmit_dir.z, 2);
        float B = 2 * ( transmit_dir.x*(start.x - spheres[j].center.x) + transmit_dir.y*(start.y - spheres[j].center.y) + transmit_dir.z*(start.z - spheres[j].center.z) );
        float C = pow( start.x - spheres[j].center.x, 2) + pow( start.y - spheres[j].center.y, 2) + pow( start.z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
        float D = B*B - 4*A*C;
        
        if( D < 0 ) continue;
        
        float t0 = (-B - sqrt(D))/2*A;
        float t1 = (-B + sqrt(D))/2*A;
        if( t0 < t1 ) local_t = t0; else local_t = t1;
        
        if( local_t < 0 ) continue;
        
        // farther from object
        if( local_t > t ) continue;
        
        t = local_t;
        
        ray = start + *direction * t;
        
        vector3 *cen = new vector3( spheres[j].center.x, spheres[j].center.y, spheres[j].center.z);
        vec3_unit *surface_normal = new vec3_unit(*cen, ray);
        
        //표면에서 eye로 가는 방향
        vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
        
        is_rec_calculating = true;
        
        if( in == REFLECTION )
        {
            rec_refl_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refl_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refl_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            
            is_rec_reflected = true;
        }
        else
        {
            rec_refr_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refr_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            rec_refr_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *surface_normal );
            
            is_rec_refracted = true;
        }
        which_recursive_number = j;
        which_recursive_shape = SPHERE;
    }
    
    for( int j = 0; j < no_of_triangles; j++ )
    {
        if( s == TRIANGLE && i == j ) continue;
        
        float local_t;
        
        vec3_unit pn = triangles[j].surface_normal;
        float eqa = triangles[j].eqa;
        float eqb = triangles[j].eqb;
        float eqc = triangles[j].eqc;
        float eqd = triangles[j].eqd;
        
        if( eqa*transmit_dir.x + eqb*transmit_dir.y + eqc*transmit_dir.z == 0 ) continue;
        
        local_t =  eqa*start.x + eqb*start.y + eqc*start.z + eqd;
        local_t /= eqa*transmit_dir.x + eqb*transmit_dir.y + eqc*transmit_dir.z;
        local_t = -local_t;
        
        if( local_t < 0 ) continue;
        
        // farther from object
        if( local_t > t ) continue;
        
        vec3_unit surface_normal = pn;
        
        ray = start + *direction * local_t;
        
        // containment test
        vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
        
        vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
        vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
        
        if( cp1*cp2 < 0 ) continue;
        
        vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
        vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
        
        if( cp3*cp4 < 0 ) continue;
        
        vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
        vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
        
        if( cp5*cp6 < 0 ) continue;
        
        t = local_t;
        
        // 표면에서 eye로 가는 방향
        vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
        
        is_rec_calculating = true;
        
        if( in == REFLECTION )
        {
            rec_refl_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refl_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refl_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            
            is_rec_reflected = true;
        }
        else
        {
            rec_refr_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refr_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            rec_refr_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray,  surface_normal, *view_dir,  surface_normal );
            
            is_rec_refracted = true;
        }
        
        which_recursive_number = j;
        which_recursive_shape = TRIANGLE;
    }
    
    return true;
}

void init(void)
{
    glEnable(GL_DEPTH_TEST);
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel (GL_FLAT);
    glPointSize( point_size );
    glLineWidth( 25.0 );
}

void display(void)
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(eyex, eyey, eyez, lookAt_x, lookAt_y, lookAt_z, up_x, up_y, up_z);
    
    //r += 10;
    //glRotatef(r,0,10,0);
    
    // 큰 점으로 안그려지길래 작은 원으로..
    for( int i = 0; i < no_of_light_source; i++ )
    {
        glPushMatrix();
        glTranslatef( lights[i].position.x, lights[i].position.y, lights[i].position.z);
        glColor3f( 1.0, 1.0, 1.0);
        glutSolidSphere( 0.3, 10, 10 );
        glPopMatrix();
    }
    
    // DRAW SPHERES
    for( int i = 0; i < no_of_spheres; i++ )
    {
        float tx , ty, tz;
        float theta, pi;
        
        // 중심, 반지름 저장해놓은데서 꺼내오기
        vector3 *cen = new vector3( spheres[i].center.x, spheres[i].center.y, spheres[i].center.z);
        float radius = spheres[i].radius;
        
        glPushMatrix();
        glBegin(GL_POINTS);
        
        for( theta = 0; theta < 180; theta+= theta_inc )
        {
            for( pi = 0; pi < 360; pi+= PI_inc )
            {
                // 원통형 좌표계로 원 색칠
                tx =  cen->x + radius * sin(theta) * cos( pi );
                ty =  cen->y + radius * sin(theta) * sin( pi );
                tz =  cen->z + radius * cos(theta);
                
                // 법선벡터
                vector3 *surface_position = new vector3( tx, ty, tz);
                
                // 다른 물체가 앞에 있어서 안보이는 경우
                if( !is_visible( *surface_position, SPHERE, i) ) continue;
                
                vec3_unit *surface_normal = new vec3_unit(*cen, *surface_position);
                
                vector3 *eye_position = new vector3( eyex, eyey, eyez);
                vec3_unit *view_dir = new vec3_unit( *surface_position, *eye_position);
                
                // 내 눈에서 안보이는 경우
                if( dot_product( *view_dir, *surface_normal) < 0 ) continue;
                
                // initially, no vector3 is illuminated either by diffuse or specular intensity at ( tx, ty, tz)
                // initially, no vector3 is reflecting any other object
                diffuse         = 0.0;
                specular        = 0.0;
                refl_R   = 0.0;
                refl_G = 0.0;
                refl_B  = 0.0;
                refr_R   = 0.0;
                refr_G = 0.0;
                refr_B  = 0.0;
                rec_refl_R   = 0.0;
                rec_refl_G = 0.0;
                rec_refl_B  = 0.0;
                rec_refr_R   = 0.0;
                rec_refr_G = 0.0;
                rec_refr_B  = 0.0;
                
                
                is_reflected = false;
                is_refracted = false;
                is_rec_reflected = false;
                is_rec_refracted = false;
                is_rec_calculating = false;
                
                in_shadow    = false;
                
                // LOCAL INTENSITY DUE TO LIGHT SOURCES
                
                // calculate intensities contributed by each light source
                for( int j = 0; j < no_of_light_source; j++ )
                {
                    
                    // diffuse and specular contributed by this particular light source
                    float local_diffuse = 0.0, local_specular = 0.0;
                    
                    // light source position
                    vector3 *light_source_position = new vector3( lights[j].position.x, lights[j].position.y, lights[j].position.z);
                    
                    // light direction
                    vec3_unit *light_dir = new vec3_unit( *light_source_position, *surface_position);
                    -(*light_dir);
                    
                    // SHADOW, if in shadow, check for the next light source
                    
                    if( is_shadow(*surface_position, *light_source_position, SPHERE, i) ) continue;
                    
                    // DIFFUSE
                    
                    float diffuse_dot_prod = dot_product( *light_dir, *surface_normal);
                    
                    // if (tx,ty,tz) is not visible by the light source, no need to proceed for this light source
                    if( diffuse_dot_prod < 0 ) continue;
                    
                    local_diffuse = spheres[i].kd * lights[j].intensity * diffuse_dot_prod;
                    
                    // SPECULAR
                    
                    // compute Half Way Vector H = (L+V)/|L+V|, NOTE half_way_vector is not a pointer
                    vec3_unit half_way_vector = *light_dir + *view_dir;
                    half_way_vector.unit_vectorize();
                    
                    float specular_dot_prod = pow ( dot_product( *surface_normal, half_way_vector), spheres[i].n);
                    local_specular += spheres[i].ks * lights[j].intensity * specular_dot_prod;
                    
                    //ATTENUATE THE INTENSITIES
                    
                    local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                    local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                    
                    diffuse += local_diffuse;
                    specular += local_specular;
                }
                
                
                // REFLECTION OF OTHER OBJECTS AT (tx, ty, tz)
                vec3_unit reflect_dir;
                vec3_unit temp = *surface_normal * (2 * dot_product( *view_dir, *surface_normal));
                reflect_dir = temp - *view_dir;
                reflect_dir.unit_vectorize();
                
                float t = 1000;
                vector3 *start = surface_position;
                vector3 *direction = new vector3( reflect_dir.x, reflect_dir.y, reflect_dir.z);
                
                // this is the direction of ray from (tx,ty,tz) in the direction of reflection
                vector3 ray = *start + *direction * t;
                
                for( int j = 0; j < no_of_spheres; j++ )
                {
                    // if negligible reflection coefficient, skip reflection calculations
                    if( spheres[i].refl_cof <= MIN ) break;
                    
                    if( i == j ) continue;
                    
                    float local_t;
                    
                    float A = pow( reflect_dir.x, 2) + pow( reflect_dir.y, 2) + pow( reflect_dir.z, 2);
                    float B = 2 * ( reflect_dir.x*(start->x - spheres[j].center.x) + reflect_dir.y*(start->y - spheres[j].center.y) + reflect_dir.z*(start->z - spheres[j].center.z) );
                    float C = pow( start->x - spheres[j].center.x, 2) + pow( start->y - spheres[j].center.y, 2) + pow( start->z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
                    float D = B*B - 4*A*C;
                    
                    if( D < 0 ) continue;
                    
                    float t0 = (-B - sqrt(D))/2*A;
                    float t1 = (-B + sqrt(D))/2*A;
                    if( t0 < t1 ) local_t = t0; else local_t = t1;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    t = local_t;
                    
                    ray = *start + *direction * t;
                    
                    vector3 *cen = new vector3( spheres[j].center.x, spheres[j].center.y, spheres[j].center.z);
                    vec3_unit *surface_normal = new vec3_unit(*cen, ray);
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refl_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                    refl_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                    refl_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                    
                    is_reflected = true;
                    
                }
                
                for( int j = 0; j < no_of_triangles; j++ )
                {
                    // if negligible reflection coefficient, skip reflection calculations
                    if( spheres[i].refl_cof <= MIN ) break;
                    
                    
                    float local_t;
                    
                    vec3_unit pn = triangles[j].surface_normal;
                    float eqa = triangles[j].eqa;
                    float eqb = triangles[j].eqb;
                    float eqc = triangles[j].eqc;
                    float eqd = triangles[j].eqd;
                    
                    if( eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z == 0 ) continue;
                    
                    local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                    local_t /= eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z;
                    local_t = -local_t;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    vec3_unit surface_normal = pn;
                    
                    ray = *start + *direction * local_t;
                    
                    // containment test
                    vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
                    
                    vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
                    vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
                    
                    if( cp1*cp2 < 0 ) continue;
                    
                    vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
                    vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
                    
                    if( cp3*cp4 < 0 ) continue;
                    
                    vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
                    vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
                    
                    if( cp5*cp6 < 0 ) continue;
                    
                    t = local_t;
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refl_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    refl_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    refl_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    
                    is_reflected = true;
                    
                }
                
                // REFLECTION PART OVER
                
                // REFRACTION OF OTHER OBJECTS AT (tx, ty, tz)
                // NOTE : THIS IS THE COMPLEX PART BECAUSE IT IS REFRACTION THROUGH SPHERE
                
                // 굴절광 찾기
                float cos_I    = dot_product( *view_dir, *surface_normal);
                float sin_I    = 1 - pow(cos_I, 2);
                float refr_index = 1/(spheres[i].refr_ind);
                float sin_T    = refr_index * sin_I;
                float cos_T    = 1 - pow(sin_T, 2);
                
                vec3_unit temp_I = *view_dir;
                -temp_I;
                
                vector3 *I = new vector3( temp_I.x         , temp_I.y         , temp_I.z          );
                vector3 *N = new vector3( surface_normal->x, surface_normal->y, surface_normal->z );
                
                // this is the refraction direction from surface to in the direction of refracted ray
                *direction = ( *I + *N * cos_I ) * refr_index - *N * cos_T;
                
                vec3_unit *refract_dir = new vec3_unit( new vector3(0.0,0.0,0.0), direction);
                
                // t is at infinity i.e. the object to be refracted is at infinity
                t = 1000;
                start = surface_position;
                direction = new vector3( refract_dir->x, refract_dir->y, refract_dir->z);
                
                // this is the direction of ray from (tx,ty,tz) in the direction of refraction i.e. going inside the sphere
                ray = *start + *direction * t;
                
                // now, we have to find the position from where this ray will come out of the sphere
                float local_t;
                
                float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
                float B = 2 * ( refract_dir->x*(start->x - spheres[i].center.x) + refract_dir->y*(start->y - spheres[i].center.y) + refract_dir->z*(start->z - spheres[i].center.z) );
                float C = pow( start->x - spheres[i].center.x, 2) + pow( start->y - spheres[i].center.y, 2) + pow( start->z - spheres[i].center.z, 2) - pow( spheres[i].radius, 2);
                float D = B*B - 4*A*C;
                
                if( D < 0 ) continue;
                
                float t0 = (-B - sqrt(D))/2*A;
                float t1 = (-B + sqrt(D))/2*A;
                if( t0 < t1 ) local_t = t1; else local_t = t0;
                
                if( local_t < 0 ) continue;
                
                t = local_t;
                
                ray = *start + *direction * t;
                
                vector3 *cen = new vector3( spheres[i].center.x, spheres[i].center.y, spheres[i].center.z);
                surface_normal = new vec3_unit(ray, *cen);
                
                // 표면에서 eye로
                view_dir = new vec3_unit( ray, *surface_position);
                
                cos_I    = dot_product( *view_dir, *surface_normal);
                sin_I    = 1 - pow(cos_I, 2);
                refr_index = spheres[i].refr_ind;
                sin_T    = refr_index * sin_I;
                cos_T    = 1 - pow(sin_T, 2);
                
                temp_I = *view_dir;
                -temp_I;
                
                I = new vector3( temp_I.x         , temp_I.y         , temp_I.z          );
                N = new vector3( surface_normal->x, surface_normal->y, surface_normal->z );
                
                // this is the refraction direction from surface to in the direction of refracted ray
                *direction = ( *I + *N * cos_I ) * refr_index - *N * cos_T;
                
                refract_dir = new vec3_unit( new vector3(0.0,0.0,0.0), direction);
                
                // t is at infinity i.e. the object to be refracted is at infinity
                t = 1000;
                start = surface_position;
                direction = new vector3( refract_dir->x, refract_dir->y, refract_dir->z);
                
                // this is the direction of ray from (tx,ty,tz) in the direction of refraction i.e. coming outside of the sphere
                ray = *start + *direction * t;
                
                for( int j = 0; j < no_of_spheres; j++ )
                {
                    // if negligible refraction coefficient, skip refraction calculations
                    if( spheres[i].refr_cof <= MIN ) break;
                    
                    
                    if( i == j ) continue;
                    
                    float local_t;
                    
                    float A = pow( refract_dir->x, 2) + pow( refract_dir->y, 2) + pow( refract_dir->z, 2);
                    float B = 2 * ( refract_dir->x*(start->x - spheres[j].center.x) + refract_dir->y*(start->y - spheres[j].center.y) + refract_dir->z*(start->z - spheres[j].center.z) );
                    float C = pow( start->x - spheres[j].center.x, 2) + pow( start->y - spheres[j].center.y, 2) + pow( start->z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
                    float D = B*B - 4*A*C;
                    
                    if( D < 0 ) continue;
                    
                    float t0 = (-B - sqrt(D))/2*A;
                    float t1 = (-B + sqrt(D))/2*A;
                    if( t0 < t1 ) local_t = t0; else local_t = t1;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    t = local_t;
                    
                    ray = *start + *direction * t;
                    
                    vector3 *cen = new vector3( spheres[j].center.x, spheres[j].center.y, spheres[j].center.z);
                    vec3_unit *surface_normal = new vec3_unit(*cen, ray);
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refr_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );
                    refr_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );
                    refr_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, *refract_dir );
                    
                    is_refracted = true;
                    
                }
                
                for( int j = 0; j < no_of_triangles; j++ )
                {
                    // if negligible refraction coefficient, skip refraction calculations
                    if( spheres[i].refr_cof <= MIN ) break;
                    
                    float local_t;
                    
                    vec3_unit pn = triangles[j].surface_normal;
                    float eqa = triangles[j].eqa;
                    float eqb = triangles[j].eqb;
                    float eqc = triangles[j].eqc;
                    float eqd = triangles[j].eqd;
                    
                    if( eqa*refract_dir->x + eqb*refract_dir->y + eqc*refract_dir->z == 0 ) continue;
                    
                    local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                    local_t /= eqa*refract_dir->x + eqb*refract_dir->y + eqc*refract_dir->z;
                    local_t = -local_t;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    vec3_unit surface_normal = pn;
                    
                    ray = *start + *direction * local_t;
                    
                    // containment test
                    vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
                    
                    vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
                    vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
                    
                    if( cp1*cp2 < 0 ) continue;
                    
                    vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
                    vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
                    
                    if( cp3*cp4 < 0 ) continue;
                    
                    vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
                    vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
                    
                    if( cp5*cp6 < 0 ) continue;
                    
                    t = local_t;
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refr_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );
                    refr_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );
                    refr_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, *refract_dir );
                    
                    is_refracted = true;
                    
                }
                
                // REFRACTION PART OVER
                
                if( is_reflected )
                {
                    if( is_refracted )
                    {
                        glColor3f((ambient * spheres[i].ka + diffuse) * spheres[i].cr * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + specular + refl_R   * ( spheres[i].refl_cof ) + refr_R   * ( spheres[i].refr_cof ),
                                  (ambient * spheres[i].ka + diffuse) * spheres[i].cg * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + specular + refl_G * ( spheres[i].refl_cof ) + refr_G * ( spheres[i].refr_cof ),
                                  (ambient * spheres[i].ka + diffuse) * spheres[i].cb * (1 - spheres[i].refl_cof - spheres[i].refr_cof ) + specular + refl_B  * ( spheres[i].refl_cof ) + refr_B  * ( spheres[i].refr_cof ));
                    }
                    else
                        glColor3f((ambient * spheres[i].ka + diffuse) * spheres[i].cr * (1 - spheres[i].refl_cof ) + specular + refl_R   * ( spheres[i].refl_cof ),
                                  (ambient * spheres[i].ka + diffuse) * spheres[i].cg * (1 - spheres[i].refl_cof ) + specular + refl_G * ( spheres[i].refl_cof ),
                                  (ambient * spheres[i].ka + diffuse) * spheres[i].cb * (1 - spheres[i].refl_cof ) + specular + refl_B  * ( spheres[i].refl_cof ));
                }
                else if( is_refracted )
                    glColor3f((ambient * spheres[i].ka + diffuse) * spheres[i].cr * (1 - spheres[i].refr_cof ) + specular + refr_R   * ( spheres[i].refr_cof ),
                              (ambient * spheres[i].ka + diffuse) * spheres[i].cg * (1 - spheres[i].refr_cof ) + specular + refr_G   * ( spheres[i].refr_cof ),
                              (ambient * spheres[i].ka + diffuse) * spheres[i].cb * (1 - spheres[i].refr_cof ) + specular + refr_B   * ( spheres[i].refr_cof ));
                else
                {
                    glPointSize( 3.0 );
                    glColor3f((ambient * spheres[i].ka + diffuse) * spheres[i].cr + specular,
                              (ambient * spheres[i].ka + diffuse) * spheres[i].cg + specular,
                              (ambient * spheres[i].ka + diffuse) * spheres[i].cb + specular);
                }
                
                glVertex3f( tx, ty, tz);
            }
        }
        
        glEnd();
        glPopMatrix();
        
    }
    
    // DRAW TRIANGLES
    for( int i = 0; i < no_of_triangles; i++ )
    {
        //std::cout << "T" << i << "\n";
        
        float tx , ty, tz;
        
        glPushMatrix();
        
        // end points of triangle
        vector3 v1 = triangles[i].v1;
        vector3 v2 = triangles[i].v2;
        vector3 v3 = triangles[i].v3;
        
        vec3_unit surface_normal = triangles[i].surface_normal;
        
        glBegin(GL_POINTS);
        
        int count = 0;
        for( float width = 0; width <= 1; width += width_inc )
        {
            float x_begin = width * (v3.x - v1.x) + v1.x ;
            float y_begin = width * (v3.y - v1.y) + v1.y ;
            float z_begin = width * (v3.z - v1.z) + v1.z ;
            
            float x_end = width * (v3.x - v2.x) + v2.x ;
            float y_end = width * (v3.y - v2.y) + v2.y ;
            float z_end = width * (v3.z - v2.z) + v2.z ;
            
            // draw line from v_begin to v_end
            for( float height = 0; height <= 1; height += height_inc )
            {
                glPointSize( point_size );
                
                // triangle vector3 at this point is tx,ty,tz
                tx = height * (x_end - x_begin) + x_begin ;
                ty = height * (y_end - y_begin) + y_begin ;
                tz = height * (z_end - z_begin) + z_begin ;
                
                // surface_position의 밝기 지정
                vector3 *surface_position = new vector3( tx, ty, tz);
                
                // 다른 물체에 가려서안보임
                if( !is_visible( *surface_position, TRIANGLE, i) ) continue;
                
                vector3 *eye_position = new vector3( eyex, eyey, eyez);
                vec3_unit *view_dir = new vec3_unit( *surface_position, *eye_position);
                
                //initialize
                diffuse         = 0.0;
                specular        = 0.0;
                refl_R   = 0.0;
                refl_G = 0.0;
                refl_B  = 0.0;
                refr_R   = 0.0;
                refr_G = 0.0;
                refr_B  = 0.0;
                rec_refl_R   = 0.0;
                rec_refl_G = 0.0;
                rec_refl_B  = 0.0;
                rec_refr_R   = 0.0;
                rec_refr_G = 0.0;
                rec_refr_B  = 0.0;
                
                
                is_reflected = false;
                is_refracted = false;
                is_rec_reflected = false;
                is_rec_refracted = false;
                is_rec_calculating = false;
                
                in_shadow    = false;
                
                // calculate intensities contributed by each light source
                for( int j = 0; j < no_of_light_source; j++ )
                {
                    
                    // diffuse and specular contributed by this particular light source
                    float local_diffuse = 0.0, local_specular = 0.0;
                    
                    // light source position
                    vector3 *light_source_position = new vector3( lights[j].position.x, lights[j].position.y, lights[j].position.z);
                    
                    // light direction
                    vec3_unit *light_dir = new vec3_unit( *light_source_position, *surface_position);
                    -(*light_dir);
                    
                    // SHADOW, if in shadow, check for the next light source
                    
                    if( is_shadow(*surface_position, *light_source_position, TRIANGLE, i) ) continue;
                    
                    // DIFFUSE
                    
                    float diffuse_dot_prod;
                    
                    if( dot_product( *view_dir, surface_normal) > 0 && dot_product( *light_dir, surface_normal) > 0 )
                    {
                        diffuse_dot_prod = dot_product( *light_dir, surface_normal);
                        local_diffuse = triangles[i].kd * lights[j].intensity * diffuse_dot_prod;
                    }
                    else if( dot_product( *view_dir, surface_normal) < 0 && dot_product( *light_dir, surface_normal) < 0 )
                    {
                        diffuse_dot_prod = -dot_product( *light_dir, surface_normal);
                        local_diffuse = triangles[i].kd * lights[j].intensity * diffuse_dot_prod;
                    }
                    else
                        continue;
                    
                    
                    vec3_unit half_way_vector = *light_dir + *view_dir;
                    half_way_vector.unit_vectorize();
                    
                    //SPECULAR
                    
                    float specular_dot_prod = pow ( dot_product( surface_normal, half_way_vector), triangles[i].n);
                    local_specular = triangles[i].ks * lights[j].intensity * specular_dot_prod;
                    
                    //attenuation
                    local_diffuse  /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                    local_specular /= constant_atten + (linear_atten)*distance(*light_source_position, *surface_position) + (quadratic_atten)*distance(*light_source_position, *surface_position)*distance(*light_source_position, *surface_position);
                    
                    diffuse += local_diffuse;
                    specular += local_specular;
                }
                
                
                // REFLECTION
                
                // this is the reflection direction
                vec3_unit reflect_dir;
                vec3_unit temp = surface_normal * (2 * dot_product( *view_dir, surface_normal));
                reflect_dir = temp - *view_dir;
                
                // t is at infinity i.e. the object to be reflected is at infinity
                float t = 1000;
                vector3 *start = surface_position;
                vector3 *direction = new vector3( reflect_dir.x, reflect_dir.y, reflect_dir.z);
                
                // this is the direction of ray from (tx,ty,tz) in the direction of reflection
                vector3 ray = *start + *direction * t;
                
                // 구 체크
                for( int j = 0; j < no_of_spheres; j++ )
                {
                    // if negligible reflection coefficient, skip reflection calculations
                    if( triangles[i].refl_cof <= MIN ) break;
                    
                    float local_t;
                    
                    float A = pow( reflect_dir.x, 2) + pow( reflect_dir.y, 2) + pow( reflect_dir.z, 2);
                    float B = 2 * ( reflect_dir.x*(start->x - spheres[j].center.x) + reflect_dir.y*(start->y - spheres[j].center.y) + reflect_dir.z*(start->z - spheres[j].center.z) );
                    float C = pow( start->x - spheres[j].center.x, 2) + pow( start->y - spheres[j].center.y, 2) + pow( start->z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
                    float D = B*B - 4*A*C;
                    
                    if( D < 0 ) continue;
                    
                    float t0 = (-B - sqrt(D))/2*A;
                    float t1 = (-B + sqrt(D))/2*A;
                    if( t0 < t1 ) local_t = t0; else local_t = t1;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    t = local_t;
                    
                    ray = *start + *direction * t;
                    
                    vector3 *cen = new vector3( spheres[j].center.x, spheres[j].center.y, spheres[j].center.z);
                    vec3_unit *surface_normal = new vec3_unit(*cen, ray);
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refl_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir );
                    refl_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir  );
                    refl_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, reflect_dir  );
                    
                    is_reflected = true;
                    
                }
                
                //삼각형 체크
                for( int j = 0; j < no_of_triangles; j++ )
                {
                    // 얼마나 돌릴건지 결정
                    if( triangles[i].refl_cof <= MIN) break;
                    
                    if( i == j ) continue;
                    
                    float local_t;
                    
                    vec3_unit pn = triangles[j].surface_normal;
                    float eqa = triangles[j].eqa;
                    float eqb = triangles[j].eqb;
                    float eqc = triangles[j].eqc;
                    float eqd = triangles[j].eqd;
                    
                    if( eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z == 0 ) continue;
                    
                    local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                    local_t /= eqa*reflect_dir.x + eqb*reflect_dir.y + eqc*reflect_dir.z;
                    local_t = -local_t;
                    
                    if( local_t < 0 ) continue;
                    
                    // farther from object
                    if( local_t > t ) continue;
                    
                    vec3_unit surface_normal = pn;
                    
                    ray = *start + *direction * local_t;
                    
                    // containment test
                    vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
                    
                    vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
                    vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
                    
                    if( cp1*cp2 < 0 ) continue;
                    
                    vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
                    vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
                    
                    if( cp3*cp4 < 0 ) continue;
                    
                    vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
                    vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
                    
                    if( cp5*cp6 < 0 ) continue;
                    
                    t = local_t;
                    
                    //surface -> eye
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refl_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    refl_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    refl_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, reflect_dir );
                    
                    is_reflected = true;
                    
                }
                
                // REFLECTION 끝
                
                // REFRACTION OF OTHER OBJECTS AT (tx, ty, tz)
                
                vec3_unit refract_dir = *view_dir;
                -refract_dir;
                
                //initialization ->일단 큰값
                t = 1000;
                start = surface_position;
                direction = new vector3( refract_dir.x, refract_dir.y, refract_dir.z);
                
                // refraction direction ray
                ray = *start + *direction * t;
                
                for( int j = 0; j < no_of_spheres; j++ )
                {
                    //control # of bounce
                    if( triangles[i].refr_cof <= MIN) break;
                    
                    float local_t;
                    
                    float A = pow( refract_dir.x, 2) + pow( refract_dir.y, 2) + pow( refract_dir.z, 2);
                    float B = 2 * ( refract_dir.x*(start->x - spheres[j].center.x) + refract_dir.y*(start->y - spheres[j].center.y) + refract_dir.z*(start->z - spheres[j].center.z) );
                    float C = pow( start->x - spheres[j].center.x, 2) + pow( start->y - spheres[j].center.y, 2) + pow( start->z - spheres[j].center.z, 2) - pow( spheres[j].radius, 2);
                    float D = B*B - 4*A*C;
                    
                    if( D < 0 ) continue;
                    
                    float t0 = (-B - sqrt(D))/2*A;
                    float t1 = (-B + sqrt(D))/2*A;
                    if( t0 < t1 ) local_t = t0; else local_t = t1;
                    
                    if( local_t < 0 ) continue;
                    if( local_t > t ) continue;
                    
                    t = local_t;
                    ray = *start + *direction * t;
                    
                    vector3 *cen = new vector3( spheres[j].center.x, spheres[j].center.y, spheres[j].center.z);
                    vec3_unit *surface_normal = new vec3_unit(*cen, ray);
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refr_R   = calculate_light_intensity( RED,   SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );
                    refr_G = calculate_light_intensity( GREEN, SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );
                    refr_B  = calculate_light_intensity( BLUE,  SPHERE, j, ray, *surface_normal, *view_dir, refract_dir );
                    
                    is_refracted = true;
                }
                
                for( int j = 0; j < no_of_triangles; j++ )
                {
                    // if negligible refraction coefficient, skip refraction calculations
                    if( triangles[i].refr_cof <= MIN ) break;
                    
                    if( i == j ) continue;
                    
                    float local_t;
                    
                    vec3_unit pn = triangles[j].surface_normal;
                    float eqa = triangles[j].eqa;
                    float eqb = triangles[j].eqb;
                    float eqc = triangles[j].eqc;
                    float eqd = triangles[j].eqd;
                    
                    if( eqa*refract_dir.x + eqb*refract_dir.y + eqc*refract_dir.z == 0 ) continue;
                    
                    local_t =  eqa*start->x + eqb*start->y + eqc*start->z + eqd;
                    local_t /= eqa*refract_dir.x + eqb*refract_dir.y + eqc*refract_dir.z;
                    local_t = -local_t;
                    
                    if( local_t < 0 ) continue;
                    if( local_t > t ) continue;
                    
                    vec3_unit surface_normal = pn;
                    
                    ray = *start + *direction * local_t;
                    
                    // containment test
                    vector3 v1 = triangles[j].v1, v2 = triangles[j].v2, v3 = triangles[j].v3;
                    
                    vector3 cp1 = cross_product_vector3( v1-v2,  v3-v2);
                    vector3 cp2 = cross_product_vector3( v1-v2, ray-v2);
                    
                    if( cp1*cp2 < 0 ) continue;
                    
                    vector3 cp3 = cross_product_vector3( v2-v3,  v1-v3);
                    vector3 cp4 = cross_product_vector3( v2-v3, ray-v3);
                    
                    if( cp3*cp4 < 0 ) continue;
                    
                    vector3 cp5 = cross_product_vector3( v3-v1,  v2-v1);
                    vector3 cp6 = cross_product_vector3( v3-v1, ray-v1);
                    
                    if( cp5*cp6 < 0 ) continue;
                    
                    t = local_t;
                    
                    // 표면에서 eye로
                    vec3_unit *view_dir = new vec3_unit( ray, *eye_position);
                    
                    refr_R   = calculate_light_intensity( RED,   TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );
                    refr_G = calculate_light_intensity( GREEN, TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );
                    refr_B  = calculate_light_intensity( BLUE,  TRIANGLE, j, ray, surface_normal, *view_dir, refract_dir );
                    
                    is_refracted = true;
                    
                }
                
                // REFRACTION 끝
                
                // DRAW
                if( is_reflected )
                {
                    if( is_refracted )
                    {
                        glColor3f((ambient * triangles[i].ka + diffuse) * triangles[i].cr * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + specular + refl_R   * ( triangles[i].refl_cof ) + refr_R   * ( triangles[i].refr_cof ),
                                  (ambient * triangles[i].ka + diffuse) * triangles[i].cg * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + specular + refl_G * ( triangles[i].refl_cof ) + refr_G * ( triangles[i].refr_cof ),
                                  (ambient * triangles[i].ka + diffuse) * triangles[i].cb * (1 - triangles[i].refl_cof - triangles[i].refr_cof ) + specular + refl_B  * ( triangles[i].refl_cof ) + refr_B  * ( triangles[i].refr_cof ));
                    }
                    else
                        glColor3f((ambient * triangles[i].ka + diffuse) * triangles[i].cr * (1 - triangles[i].refl_cof ) + specular + refl_R   * ( triangles[i].refl_cof ),
                                  (ambient * triangles[i].ka + diffuse) * triangles[i].cg * (1 - triangles[i].refl_cof ) + specular + refl_G * ( triangles[i].refl_cof ),
                                  (ambient * triangles[i].ka + diffuse) * triangles[i].cb * (1 - triangles[i].refl_cof ) + specular + refl_B  * ( triangles[i].refl_cof ));
                }
                else if( is_refracted )
                    glColor3f((ambient * triangles[i].ka + diffuse) * triangles[i].cr * (1 - triangles[i].refr_cof ) + specular + refr_R   * ( triangles[i].refr_cof ),
                              (ambient * triangles[i].ka + diffuse) * triangles[i].cg * (1 - triangles[i].refr_cof ) + specular + refr_G   * ( triangles[i].refr_cof ),
                              (ambient * triangles[i].ka + diffuse) * triangles[i].cb * (1 - triangles[i].refr_cof ) + specular + refr_B   * ( triangles[i].refr_cof ));
                else
                {
                    glColor3f((ambient * triangles[i].ka + diffuse) * triangles[i].cr + specular,
                              (ambient * triangles[i].ka + diffuse) * triangles[i].cg + specular,
                              (ambient * triangles[i].ka + diffuse) * triangles[i].cb + specular);
                }
                glVertex3f( tx, ty, tz);
            }
        }
        
        glEnd();
        glPopMatrix();
        
    }
    
    glutSwapBuffers();
}

void reshape (int w, int h)
{
    glViewport (0, 0, w, h);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glFrustum (-width/2, width/2, -height/2, height/2, near_plane, far_plane);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
    switch( key )
    {
            
        case 'q':
            eyex += 1;
            break;
        case 'w':
            eyex -= 1;
            break;
        case 'a':
            eyez += 1;
            break;
        case 's':
            eyez -= 1;
            break;
        case 'z':
            eyey += 1;
            break;
        case 'x':
            eyey -= 1;
            break;
        case 'p':
            lights[0].position.x += 1;
            break;
        case 'o':
            lights[0].position.y += 1;
        case 'l':
            lights[0].position.z -= 1;
            break;
            
            
    }
    glutPostRedisplay();
}

int main(int argc, char** argv)
{
    FILE *fp, *fobj;
    
    //obj file
    vector3 vv;
    int i;
    char s, check;
    
    //삼각형 input
    float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    float cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind;
    
    //making file read
    fp = fopen("final.txt", "r");
    
    fscanf(fp, "%f\n\n", &point_size);
    //색칠할 좌표들 최소값
    fscanf(fp, "%f %f %f %f\n\n", &width_inc, &height_inc, &theta_inc, &PI_inc);
    fscanf(fp, "%d\n\n", &no_of_spheres);
    
    for( int i = 0; i < no_of_spheres; i++ )
    {
        float radius, x, y, z;
        
        fscanf(fp, "%f\n", &radius);
        fscanf(fp, "%f %f %f\n", &x, &y, &z);
        fscanf(fp, "%f %f %f\n", &cr, &cg, &cb);
        fscanf(fp, "%f %f %f %f\n", &ka, &kd, &ks, &n);
        fscanf(fp, "%f %f %f\n\n", &refl_cof, &refr_cof, &refr_ind);
        
        vector3 *center = new vector3( x, y, z);
        sphere *temp = new sphere( radius, *center, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind);
        spheres.push_back( *temp );
    }
    
    fscanf(fp, "%d\n\n", &no_of_triangles);
    for( int i = 0; i < no_of_triangles; i++ )
    {
        //삼각형 좌표
        fscanf(fp, "%f %f %f\n", &x1, &y1, &z1);
        fscanf(fp, "%f %f %f\n", &x2, &y2, &z2);
        fscanf(fp, "%f %f %f\n", &x3, &y3, &z3);
        
        //삼각형 속성들(색깔, K값, 상수)
        fscanf(fp, "%f %f %f\n", &cr, &cg, &cb);
        fscanf(fp, "%f %f %f %f\n", &ka, &kd, &ks, &n);
        fscanf(fp, "%f %f %f\n\n", &refl_cof, &refr_cof, &refr_ind);
        
        vector3 *v1 = new vector3(x1,y1,z1), *v2 = new vector3(x2,y2,z2), *v3 = new vector3(x3,y3,z3);
        
        //normal vector찾기 -> 평면방정식 만들려고
        vec3_unit *plane_vector_1 = new vec3_unit( *v1, *v2);
        vec3_unit *plane_vector_2 = new vec3_unit( *v1, *v3);
        vec3_unit surface_normal = cross_product(*plane_vector_2, *plane_vector_1);
        surface_normal.unit_vectorize();
        
        float eqa = surface_normal.x, eqb = surface_normal.y, eqc = surface_normal.z, eqd;
        eqd = - ( eqa*x1 + eqb*y1 + eqc*z1 );
        
        triangle *temp = new triangle( eqa, eqb, eqc, eqd, surface_normal, v1, v2, v3, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind);
        
        triangles.push_back( *temp );
    }
    
    fscanf(fp, "%f %f %f\n", &eyex, &eyey, &eyez);
    fscanf(fp, "%f %f %f\n", &lookAt_x, &lookAt_y, &lookAt_z);
    fscanf(fp, "%f %f %f\n\n", &up_x, &up_y, &up_z);
    fscanf(fp, "%f\n", &near_plane);
    fscanf(fp, "%f\n", &far_plane);
    fscanf(fp, "%f %f\n\n", &width, &height);
    fscanf(fp, "%d\n", &no_of_light_source);
    
    for( int i = 0; i < no_of_light_source; i++ )
    {
        float x ,y ,z, intensity;
        
        fscanf(fp, "%f %f %f %f\n\n", &x, &y, &z, &intensity);
        vector3 *vtemp = new vector3( x, y, z );
        light *temp = new light( *vtemp, intensity );
        lights.push_back( *temp );
    }
    
    fscanf(fp, "%f %f %f\n", &constant_atten, &linear_atten, &quadratic_atten);
    fscanf(fp, "%f\n", &ambient);
    
    
/*obj file read
     fobj = fopen("dogHouse.obj", "r");
     while (!feof(fobj)){
        fscanf(fobj, "%c", &s);
        if (s == 'v')
        {
            fscanf(fobj, "%c", &check);
            if (check == ' ')
            {
                fscanf(fobj, "%f %f %f\n", &vv.x, &vv.y, &vv.z);
                v.push_back(vv);
            }
        }
        else if (s == 'f')
        {
            fscanf(fobj, " %d//%*d %d//%*d %d//%*d\n", &ff.a1, &ff.b1, &ff.c1);
            f.push_back(ff);
        }
     }
     
     for (i = 0 ; i < f.size(); i++){
        x1 = v[f[i].a1 - 1].x, y1 = v[f[i].a1 - 1].y, z1 = v[f[i].a1 - 1].z;
        x2 = v[f[i].b1 - 1].x, y2 = v[f[i].b1 - 1].y, z2 = v[f[i].b1 - 1].z;
        x3 = v[f[i].c1 - 1].x, y3 = v[f[i].c1 - 1].y, z3 = v[f[i].c1 - 1].z;
     
        cr = 1.0, cg = 1.0, cb = 1.0;
        ka = 0.2, kd = 0.5, ks = 0.7, n = 100;
        refl_cof = 0.5, refr_cof = 0; refr_ind = 1.0;
     
     
        //여기서부터 안먹히뮤
        vector3 *v1 = new vector3(x1-10, y1+2, z1-2), *v2 = new vector3(x2-10, y2+2, z2-2), *v3 = new vector3(x3-10, y3+2, z3-2);
     
        // find normal to this triangle
        vec3_unit *plane_vector_1 = new vec3_unit( *v1, *v2);
        vec3_unit *plane_vector_2 = new vec3_unit( *v1, *v3);
        vec3_unit surface_normal = cross_product(*plane_vector_2, *plane_vector_1);
        surface_normal.unit_vectorize();
     
        float eqa = surface_normal.x, eqb = surface_normal.y, eqc = surface_normal.z, eqd;
        eqd = - ( eqa*x1 + eqb*y1 + e[i]qc*z1 );
     
        triangle *temp = new triangle( eqa, eqb, eqc, eqd, surface_normal, v1, v2, v3, cr, cg, cb, ka, kd, ks, n, refl_cof, refr_cof, refr_ind );
        triangles.push_back( *temp );
    }
*/
    
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DEPTH | GLUT_RGB);
    glutInitWindowSize (800, 450);
    glutInitWindowPosition (100, 100);
    glutCreateWindow (argv[0]);
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}
