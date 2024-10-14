//-------------------------------------------------------------------------------------------
// File : main.cpp
// Desc : expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)
//        Original Code by T.Hachisuka (https://cs.uwaterloo.ca/~thachisu/smallppm_exp.cpp)
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
// Includes
//-------------------------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <bitmap.h>
#include <math_util.h>
#include <list>
#include <vector>
#include <chrono>
#include <sphere.h>
#include <hitrecord.h>

#define PHOTON_COUNT_MUTIPLIER  1000

namespace /* anonymous */ {

//-------------------------------------------------------------------------------------------
// Constant Values
//-------------------------------------------------------------------------------------------
static const double ALPHA = 2.0 / 3.0; // the alpha parameter of PPM


//-------------------------------------------------------------------------------------------
// Global Variables
//-------------------------------------------------------------------------------------------
std::list<HitRecord*>               hitpoints;
std::vector<std::list<HitRecord*> > hash_grid;
double                              hash_s;
BoundingBox                         hpbbox;
SphereObject sph[] = 
{
    // Scene: radius, position, color, material
    SphereObject(1e5,  Vector3( 1e5 + 1,   40.8,       81.6 ), Vector3( 0.25, 0.75, 0.25 ),  MaterialType::Matte ),   //Right
    SphereObject(1e5,  Vector3(-1e5 + 99,  40.8,       81.6 ), Vector3( 0.25, 0.25, 0.75 ),  MaterialType::Matte ),   //Left
    SphereObject(1e5,  Vector3( 50,        40.8,        1e5 ), Vector3( 0.75, 0.75, 0.75 ),  MaterialType::Matte ),   //Back
    SphereObject(1e5,  Vector3( 50,        40.8, -1e5 + 170 ), Vector3( 0.01, 0.01, 0.01 ),  MaterialType::Matte ),   //Front
    SphereObject(1e5,  Vector3( 50,         1e5,       81.6 ), Vector3( 0.75, 0.75, 0.75 ),  MaterialType::Matte ),   //Bottomm
    SphereObject(1e5,  Vector3( 50, -1e5 + 81.6,       81.6 ), Vector3( 0.75, 0.75, 0.75 ),  MaterialType::Matte ),   //Top
    SphereObject(16.5, Vector3( 27,        16.5,         47 ), Vector3( 0.75, 0.25, 0.25 ),  MaterialType::Mirror),   //Mirror
    SphereObject(16.5, Vector3( 73,        16.5,         88 ), Vector3( 0.99, 0.99, 0.99 ),  MaterialType::Glass ),   //Glass
    SphereObject(8.5,  Vector3( 50,         8.5,         60 ), Vector3( 0.75, 0.75, 0.75 ),  MaterialType::Matte ),   //Middle
};

} // namespace /* anonymous */


//-------------------------------------------------------------------------------------------
//      ハッシュキーを取得します.
//-------------------------------------------------------------------------------------------
inline unsigned int hash
(
    const int ix,
    const int iy,
    const int iz
) 
{
    return (unsigned int)(
        (ix * 73856093) ^
        (iy * 19349663) ^
        (iz * 83492791)) % hash_grid.size();
}

//-------------------------------------------------------------------------------------------
//      ハッシュグリッドを構築します.
//-------------------------------------------------------------------------------------------
void build_hash_grid
(
    const int w,
    const int h
)
{
    // heuristic for initial radius
    auto size = hpbbox.maxi - hpbbox.mini;
    auto irad = ((size.x + size.y + size.z) / 3.0) / ((w + h) / 2.0) * 2.0;

    // determine hash table size
    // we now find the bounding box of all the measurement points inflated by the initial radius
    hpbbox.reset();
    auto photon_count = 0;
    for( auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr )
    {
        auto hp  = (*itr);
        hp->r2   = irad *irad;
        hp->n    = 0;
        hp->flux = Vector3();

        photon_count++;
        hpbbox.merge( hp->pos - irad );
        hpbbox.merge( hp->pos + irad );
    }

    // make each grid cell two times larger than the initial radius
    hash_s = 1.0 / (irad * 2.0);

    // build the hash table
    hash_grid.resize( photon_count );
    hash_grid.shrink_to_fit();
    for( auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr )
    {
        auto hp = (*itr);
        auto min = ((hp->pos - irad) - hpbbox.mini) * hash_s;
        auto max = ((hp->pos + irad) - hpbbox.mini) * hash_s;

        auto minX = abs(int(min.x));
        auto maxX = abs(int(max.x));

        auto minY = abs(int(min.y));
        auto maxY = abs(int(max.y));

        auto minZ = abs(int(min.z));
        auto maxZ = abs(int(max.z));

        for (int iz = minZ; iz <= maxZ; iz++)
        {
            for (int iy = minY; iy <= maxY; iy++)
            {
                for (int ix = minX; ix <= maxX; ix++)
                {
                    int hv = hash( ix, iy, iz );
                    hash_grid[ hv ].push_back( hp );
                }
            }
        }
    }
}

//-------------------------------------------------------------------------------------------
//      交差判定を行います.
//-------------------------------------------------------------------------------------------
inline bool intersect( const Ray &r, double &t, int &id )
{
    int n = sizeof(sph) / sizeof(sph[0]);
    auto d = D_INF;
    t = D_INF;
    for (int i = 0; i < n; i++)
    {
        d = sph[i].intersect(r);
        if (d < t)
        {
            t = d;
            id = i;
        }
    }

    return ( t < D_INF );
}

//-------------------------------------------------------------------------------------------
//      フォトンレイを生成します.
//-------------------------------------------------------------------------------------------
void genp( Ray* pr, Vector3* f, int i )
{
    // generate a photon ray from the point light source with QMC
    (*f) = Vector3( 2500, 2500, 2500 ) * ( D_PI * 4.0 ); // flux
    auto p  = 2.0 * D_PI * halton( 0, i );
    auto t  = 2.0 * acos( sqrt(1. - halton( 1, i ) ));
    auto st = sin( t );

    pr->dir = Vector3( cos( p ) * st, cos( t ), sin( p ) * st );
    pr->pos = Vector3( 50, 60, 85 );
}

//-------------------------------------------------------------------------------------------
//      レイを追跡します.
//-------------------------------------------------------------------------------------------
void trace( const Ray &r, int dpt, bool m, const Vector3 &fl, const Vector3 &adj, int i )
{
    double t;
    int id;

    dpt++;
    if (!intersect(r, t, id) || (dpt >= 20))
        return;

    auto d3 = dpt * 3;
    const auto &obj = sph[ id ];
    auto x  = r.pos + r.dir*t, n = normalize( x - obj.pos );
    auto f  = obj.color;
    auto nl = ( dot(n, r.dir ) < 0 ) ? n : n*-1;
    auto p  = ( f.x > f.y && f.x > f.z ) ? f.x : ( f.y > f.z ) ? f.y : f.z;

    if ( obj.type == MaterialType::Matte )
    {
        if (m) 
        {
            // eye ray
            // store the measurment point
            auto hp = new HitRecord;
            hp->f = mul(f,adj);
            hp->pos = x;
            hp->nrm = n;
            hp->idx = i;
            hitpoints.push_back( hp );

            // find the bounding box of all the measurement points
            hpbbox.merge( x );
        }
        else
        {
            // photon ray
            // find neighboring measurement points and accumulate flux via progressive density estimation
            auto hh = (x - hpbbox.mini) * hash_s;
            auto ix = abs(int(hh.x));
            auto iy = abs(int(hh.y));
            auto iz = abs(int(hh.z));
            // strictly speaking, we should use #pragma omp critical here.
            // it usually works without an artifact due to the fact that photons are 
            // rarely accumulated to the same measurement points at the same time (especially with QMC).
            // it is also significantly faster.
            {
                auto& list = hash_grid[ hash( ix, iy, iz ) ];
                for( auto itr = list.begin(); itr != list.end(); itr++ )
                {
                    auto hp = (*itr);
                    auto v = hp->pos - x;
                    // check normals to be closer than 90 degree (avoids some edge brightning)
                    if ((dot(hp->nrm,n) > 1e-3) && (dot(v,v) <= hp->r2))
                    {
                        // unlike N in the paper, hp->n stores "N / ALPHA" to make it an integer value
                        auto g = (hp->n * ALPHA + ALPHA ) / ( hp->n * ALPHA + 1.0 );
                        hp->r2 = hp->r2 * g;
                        hp->n++;
                        hp->flux = ( hp->flux + mul( hp->f, fl ) / D_PI ) * g;
                    }
                }
            }

            // use QMC to sample the next direction
            auto r1  = 2.0 * D_PI * halton( d3 - 1, i );
            auto r2  = halton( d3 + 0, i );
            auto r2s = sqrt( r2 );
            auto w   = nl;
            auto u   = normalize(cross((fabs(w.x) > .1 ? Vector3(0, 1, 0) : Vector3(1, 0, 0)), w));
            auto v   = cross( w, u );
            auto d   = normalize( u * cos( r1 ) * r2s + v * sin( r1 ) * r2s + w * sqrt( 1 - r2 ));

            if ( halton( d3 + 1, i ) < p )
                trace(Ray(x, d), dpt, m, mul(f,fl)*(1. / p), mul(f, adj), i);
        }

    }
    else if ( obj.type == MaterialType::Mirror )
    {
        trace(Ray(x, reflect(r.dir, n)), dpt, m, mul(f,fl), mul(f,adj), i);
    }
    else 
    {
        Ray lr( x, reflect( r.dir, n ) );
        auto into  = dot(n, nl ) > 0.0;
        auto nc    = 1.0;
        auto nt    = 1.5;
        auto nnt   = (into) ? nc / nt : nt / nc;
        auto ddn   = dot( r.dir, nl );
        auto cos2t = 1 - nnt * nnt * ( 1 - ddn * ddn );

        // total internal reflection
        if (cos2t < 0)
            return trace(lr, dpt, m, mul(f, fl), mul(f, adj), i);

        auto td = normalize(r.dir * nnt - n * ( ( into ? 1 : -1 ) * ( ddn * nnt + sqrt( cos2t ))));
        auto a  = nt - nc;
        auto b  = nt + nc;
        auto R0 = a * a / ( b * b );
        auto c  = 1 - (into ? -ddn : dot(td, n));
        auto Re = R0 + (1 - R0) * c * c * c * c * c;
        auto P  = Re;
        Ray  rr(x, td);
        auto fa  = mul( f, adj );
        auto ffl = mul( f, fl  );

        if (m) 
        {
            // eye ray (trace both rays)
            trace( lr, dpt, m, ffl, fa * Re, i );
            trace( rr, dpt, m, ffl, fa * (1.0 - Re), i );
        }
        else 
        {
            // photon ray (pick one via Russian roulette)
            ( halton( d3 - 1, i ) < P ) 
                ? trace( lr, dpt, m, ffl, fa * Re, i )
                : trace( rr, dpt, m, ffl, fa * (1.0 - Re), i );
        }
    }
}

//-------------------------------------------------------------------------------------------
//      eye rayを追跡します.
//-------------------------------------------------------------------------------------------
void trace_ray( int w, int h )
{
    // 開始時刻.
    auto start = std::chrono::system_clock::now();

    // trace eye rays and store measurement points
    Ray cam(
        Vector3(50, 48, 295.6),
        normalize(Vector3(0, -0.042612, -1))
    );
    auto cx = Vector3( w * 0.5135 / h, 0, 0 );
    auto cy = normalize( cross( cx, cam.dir ) ) * 0.5135;

    for (int y = 0; y < h; y++)
    {
        // 進捗率を表示.
        fprintf( stdout, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1) );

        for (int x = 0; x < w; x++) 
        {
            auto idx = x + y * w;
            auto d   = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.dir;
            trace( Ray(cam.pos + d * 140, normalize(d)), 0, true, Vector3(), Vector3(1, 1, 1), idx );
        }
    }

    // 改行.
    fprintf( stdout, "\n" );

    // 処理時間を表示.
    auto end = std::chrono::system_clock::now();
    auto dif = end - start;
    fprintf( stdout, "Ray Tracing Pass : %lld(msec)\n", std::chrono::duration_cast<std::chrono::milliseconds>(dif).count() );

    start = std::chrono::system_clock::now();

    // build the hash table over the measurement points
    build_hash_grid( w, h );

    // 処理時間を表示.
    end = std::chrono::system_clock::now();
    dif = end - start;
    fprintf( stdout, "Build Hash Grid  : %lld(msec)\n", std::chrono::duration_cast<std::chrono::milliseconds>(dif).count() );
}

//-------------------------------------------------------------------------------------------
//      photon rayを追跡します.
//-------------------------------------------------------------------------------------------
void trace_photon( int s )
{
    // 開始時刻.
    auto start = std::chrono::system_clock::now();

    auto vw = Vector3(1, 1, 1);

    // trace photon rays with multi-threading
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < s; i++) 
    {
        // 進捗率表示.
        auto p = 100.0 * ( i + 1 ) / s;
        fprintf( stdout, "\rPhotonPass %5.2f%%", p );

        int m = PHOTON_COUNT_MUTIPLIER * i;
        Ray r;
        Vector3 f;
        for ( int j = 0; j < PHOTON_COUNT_MUTIPLIER; j++ )
        {
            genp( &r, &f, m + j );
            trace( r, 0, false, f, vw, m + j );
        }
    }

    // 改行.
    fprintf( stdout, "\n" );

    // 処理時間を表示.
    auto end = std::chrono::system_clock::now();
    auto dif = end - start;
    fprintf( stdout, "Photon Tracing Pass : %lld(sec)\n", std::chrono::duration_cast<std::chrono::seconds>(dif).count() );
}

//-------------------------------------------------------------------------------------------
//      密度推定を行います.
//-------------------------------------------------------------------------------------------
void density_estimation( Vector3* color, int num_photon )
{
    // density estimation
    for( auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr )
    {
        auto hp = (*itr);
        auto i = hp->idx;
        color[i] = color[i] + hp->flux * ( 1.0 / ( D_PI * hp->r2 * num_photon * PHOTON_COUNT_MUTIPLIER ));
    }
}

//-------------------------------------------------------------------------------------------
//      メインエントリーポイントです.
//-------------------------------------------------------------------------------------------
int main(int argc, char **argv) 
{
    auto w = 1280;      // 画像の横幅.
    auto h = 1080;      // 画像の縦幅.
    auto s = 10000;     // s * 1000 個のフォトンがトレーシングされることになる (s * PHOTON_COUNT_MULTIPLIER).
    auto c = new Vector3[ w * h ];

    hpbbox.reset();

    // カメラトレーシング.
    trace_ray( w, h );

    // フォトントレーシング.
    trace_photon( s );

    // 密度推定・放射輝度推定.
    density_estimation( c, s );

    // 画像に出力.
    const auto kGamma = 2.2;
    save_to_bmp( "image.bmp", w, h, &c[0].x, kGamma );

    delete [] c;
    c = nullptr;

    return 0;
}
