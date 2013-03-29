//#define DEBUG
#include <stdio.h> 
#include <osg/Notify>
#include <osg/Geometry>
#include <osg/Geode>
#include <osg/StateSet>
#include <osg/PolygonMode>
#include <osg/Material>

#include <osgDB/WriteFile>

#include <osgUtil/DelaunayTriangulator>
#include <osgUtil/IntersectVisitor>
#include "TerrainProcessingUtils.h"
#include "CorrectionSpaceMesh.h"
#include "DelaunayWAR.h"
#include "bucket.h"

extern osg::ref_ptr<osg::Group> rootNode;		//求高度用的全局group，也可以用来求交点。
extern int LoadRootNodeByVertex(osg::Vec3 v);
extern double GetElevation(double lon, double lat, double pre);
extern double GetElevationFromHgt(double lon, double lat);
extern char tempDirectory[_MAX_PATH];
extern char debugFile[_MAX_PATH];
LimitValue currMes;

#define sqr(x)  ((x)*(x))

/* 根据一个点，转换成经纬度坐标系，求出最大最小点 */
int FindTheMaxMinValues(osg::Vec3 v)
{
	Carton ct;	ct.x = v[0]; ct.y = v[1]; ct.z = v[2];
	Geodesy gd;
	CartToGeod(&ct, &gd);
	
	/* 求取最大、最小经纬度值 */
	if(gd._lon > currMes.maxLongi) currMes.maxLongi = gd._lon;
	if(gd._lon < currMes.minLongi) currMes.minLongi = gd._lon;
	if(gd._lat > currMes.maxLati ) currMes.maxLati  = gd._lat;
	if(gd._lat < currMes.minLati ) currMes.minLati  = gd._lat;

	return 1;
}

CorrectionSpaceMesh::CorrectionSpaceMesh( osg::Group *targets, double K ):
    _valid(false),
    _K(K)
{
	osg::ref_ptr<osg::Vec3Array> csmCoords = new osg::Vec3Array;
	currMes.maxLongi = currMes.maxLati = -FLT_MAX; currMes.minLongi = currMes.minLati = FLT_MAX; 
	
	/* targets 组内实际上只有一个geode，这个geode就是机场边缘的线环。*/
	/* 求出线环上各个点与地景面垂直距离的最大值, 存在d里面 */
	double d =  0.0;
	for( unsigned int i = 0; i < targets->getNumChildren(); i++ )
	{
		osg::Geode *gnode = dynamic_cast<osg::Geode*>(targets->getChild(i));
		osg::Geode::DrawableList dl = gnode->getDrawableList();
		for(osg::Geode::DrawableList::iterator iter = dl.begin(); iter != dl.end(); iter++)
		{
			osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry*> ((*iter).get());
			osg::Vec3Array *bh = dynamic_cast<osg::Vec3Array*>(gm->getVertexArray());
			for(osg::Vec3Array::iterator p = bh->begin(); p != bh->end(); p++)
			{
				/* 根据*p来求出经纬度的最大、最小值。 */
				FindTheMaxMinValues(*p);
				
				/* 计算每个点同地景面的垂直相交点，并求出距离这个相交点的最大值 */
				osg::Vec3 isect;
				while( TerrainProcessingUtils::getIntersect( rootNode, *p, isect ) == false )
				{
					/* 根据*p来确定rootNode，并同时求出经纬度的最大、最小值。 */
					if(!LoadRootNodeByVertex(*p))
					{
						/* 换一种方式求交点 */
						Geodesy geod;
						Carton  cart; cart.x = p->x(); cart.y = p->y(); cart.z = p->z();
						CartToGeod(&cart, &geod);
						double elev = GetElevationFromHgt(geod._lon, geod._lat);
						if(elev == 1.0)
						{
							osg::notify(osg::WARN) << "Terrain Deformation Software:  ERROR in Target Database Base definition:  at least one target is not completely over terrain." << std::endl;
							return;
						}
						geod._elevation = elev;
						GeodToCart(&geod, &cart);
						isect.set(cart.x, cart.y, cart.z);
						break;
					}
				}
				
				double l = ((*p) - isect).length();
				if( l > d ) d = l;
				
				csmCoords->push_back( *p );
			}
		}
	}

	/* 求出四个极值点的经纬高，然后算出之间的距离 */
	Geodesy G00, G01, G10, G11;
	G00._lon = currMes.minLongi; G00._lat = currMes.minLati; G00._elevation = GetElevation(G00._lon, G00._lat, 0.0);
	G01._lon = currMes.minLongi; G01._lat = currMes.maxLati; G01._elevation = GetElevation(G01._lon, G01._lat, 0.0);
	G10._lon = currMes.maxLongi; G10._lat = currMes.minLati; G10._elevation = GetElevation(G10._lon, G10._lat, 0.0);
	G11._lon = currMes.maxLongi; G11._lat = currMes.maxLati; G11._elevation = GetElevation(G11._lon, G11._lat, 0.0);
	double distLonM, distLatM;
	distLonM = distanceM(&G00, &G10);	distLatM = distanceM(&G00, &G01);

	/* 最大垂直距离放大固定的倍数，作为水平扩展距离 */
	d *= _K; // Correction should be zero if more than (d * _K) away from any target vertex.
	double extLon, extLat;
	extLon = (d / distLonM) * (currMes.maxLongi - currMes.minLongi);
	extLat = (d / distLatM) * (currMes.maxLati  - currMes.minLati );
	
	currMes.maxLongi += extLon;
	currMes.minLongi -= extLon;
	currMes.maxLati  += extLat;
	currMes.minLati  -= extLat;

	G00._lon = currMes.minLongi; G00._lat = currMes.minLati; G00._elevation = GetElevation(G00._lon, G00._lat, 0.0);
	G01._lon = currMes.minLongi; G01._lat = currMes.maxLati; G01._elevation = GetElevation(G01._lon, G01._lat, 0.0);
	G10._lon = currMes.maxLongi; G10._lat = currMes.minLati; G10._elevation = GetElevation(G10._lon, G10._lat, 0.0);
	G11._lon = currMes.maxLongi; G11._lat = currMes.maxLati; G11._elevation = GetElevation(G11._lon, G11._lat, 0.0);

	Carton C00, C01, C10, C11;
	osg::ref_ptr<osg::Vec3Array >periph = new osg::Vec3Array;
	GeodToCart(&G00, &C00);		periph->push_back( osg::Vec3( C00.x, C00.y, C00.z ));
	GeodToCart(&G01, &C01);		periph->push_back( osg::Vec3( C01.x, C01.y, C01.z ));
	GeodToCart(&G10, &C10);		periph->push_back( osg::Vec3( C10.x, C10.y, C10.z ));
	GeodToCart(&G11, &C11);		periph->push_back( osg::Vec3( C11.x, C11.y, C11.z ));

	for( osg::Vec3Array::iterator p = periph->begin(); p != periph->end(); p++ )
		csmCoords->push_back( *p );

	///////////////////////////////////////////////////////
	// If needed only ....
	osg::ref_ptr<DelaunayWAR> dt = new DelaunayWAR(csmCoords.get() );
	//osg::ref_ptr<osgUtil::DelaunayTriangulator> dt = new osgUtil::DelaunayTriangulator(csmCoords.get() );
	dt->triangulate();
	
	osg::DrawElementsUInt *t = dt->getTriangles();
	osg::DrawElementsUInt::iterator p;
	
	for( p = t->begin() ; p != t->end(); p+=3 )
	{
		Triangle tri((*csmCoords)[*(p+0)], (*csmCoords)[*(p+1)], (*csmCoords)[*(p+2)]); 
		tri.setK( _K );
		_triangles.push_back( tri );
	}

	_valid = true;
	
	// Dump out CSM as csm.osg for debug purposes.
	calcCSMGeode();
}

void CorrectionSpaceMesh::calcCSMGeode()
{
	if( !_valid )   return;
	
	osg::ref_ptr<osg::Vec3Array>coords = new osg::Vec3Array;
	for( std::vector<Triangle>::iterator p =  _triangles.begin();
	        p != _triangles.end(); p++ )
	{
		osg::Vec3 A, B, C;
		p->getABC( A, B, C );
		coords->push_back( A );
		coords->push_back( B );
		coords->push_back( C );
	}
	
	osg::ref_ptr<osg::Geometry> g = new osg::Geometry;
	g->setVertexArray( coords.get() );
	g->addPrimitiveSet( new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES, 0, coords->size()));
	
	osg::PolygonMode *pm = new osg::PolygonMode( osg::PolygonMode::FRONT_AND_BACK, osg::PolygonMode::LINE );
	osg::StateSet *sset = new osg::StateSet;
	sset->setAttributeAndModes( pm );
	sset->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
	g->setStateSet( sset );
	
	csmGd = new osg::Geode;
	csmGd->addDrawable( g.get() );
	
	/* 计算包围球 */
	osg::BoundingSphere bbx;
	bbx = csmGd->getBound();
	center = bbx.center();
	radius = bbx.radius();
}

osg::Vec3 CorrectionSpaceMesh::getZCorrection( const osg::Vec3 &v )
{
	if( !_valid )		return v;

	Geodesy gd;
	Carton  ct;	
	osg::Vec3d gv;
	ct.x = v[0]; ct.y = v[1]; ct.z = v[2];	CartToGeod(&ct, &gd); 	gv.set(gd._lon, gd._lat, gd._elevation);

	std::vector<Triangle>::iterator p;
	for( p = _triangles.begin(); p != _triangles.end(); p++ )
	{
		if( p->within( gv ) )
		{
			gd._elevation = p->C(gv);
			GeodToCart(&gd, &ct);
			osg::Vec3 v; v.set(ct.x, ct.y, ct.z);
			return v;
		}
	}

	return v;
}


void CorrectionSpaceMesh::correctPoint( osg::Vec3 &p, bool print )
{
	osg::Vec3 z = getZCorrection(p);
	p = z;
}

CorrectionSpaceMesh::Triangle::Triangle( osg::Vec3 _a, osg::Vec3 _b, osg::Vec3 _c ):
    a(_a),
    b(_b),
    c(_c)
{
	/* 根据直角坐标系的值来计算经纬高坐标系的值 */
	Geodesy gd;
	Carton  ct;
	
	ct.x = a[0]; ct.y = a[1]; ct.z = a[2];	CartToGeod(&ct, &gd); 	x.set(gd._lon, gd._lat, gd._elevation);
	ct.x = b[0]; ct.y = b[1]; ct.z = b[2];	CartToGeod(&ct, &gd); 	y.set(gd._lon, gd._lat, gd._elevation);
	ct.x = c[0]; ct.y = c[1]; ct.z = c[2];	CartToGeod(&ct, &gd); 	z.set(gd._lon, gd._lat, gd._elevation);
}

void CorrectionSpaceMesh::Triangle::getABC( osg::Vec3 &A, osg::Vec3 &B, osg::Vec3 &C )
{
	A = a;
	B = b;
	C = c;
}

void CorrectionSpaceMesh::Triangle::print()
{
	printf(" %f %f %f  -  %f %f %f  -  %f %f %f\n",
		a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
}

bool CorrectionSpaceMesh::Triangle::within( double m, double n )
{
	osg::Vec3d p(m,n,0.0);
	
	if( _det( x, y, p ) < 0 )	return false;
	if( _det( y, z, p ) < 0 )	return false;
	if( _det( z, x, p ) < 0 )	return false;
	return true;
}

bool CorrectionSpaceMesh::Triangle::within( const osg::Vec3d &v )
{
	return within( v[0], v[1] );
}

double CorrectionSpaceMesh::Triangle::_det( osg::Vec3d A, osg::Vec3d B, osg::Vec3d P )
{
	return (B[0] - A[0]) * (P[1] - A[1]) - (P[0] - A[0]) * (B[1] - A[1]);
}

double CorrectionSpaceMesh::Triangle::f(osg::Vec3d p, double m, double n )
{
	double ai = p[2];
	double d_2 = sqr(m - p[0]) + sqr(n - p[1]);
	
	double ff;
	if( ai != 0.0 )
		ff = -d_2/(sqr(_K)*sqr(ai));
	else
		return 0;

	return ai * exp(ff);
}


// Height Correction Function  
// C(x,y) = w1*f1(x,y) + w2*f2(x,y) + w3*f3(x,y)
double CorrectionSpaceMesh::Triangle::C( const osg::Vec3d &p )
{
	osg::Vec3d I;
	double w1, w2, w3;
	double l;
	
	_intersect(x, p, y, z, I );
	l = _length2D(x,I);
	if( l != 0.0 )
		w1 = 1.0 - (_length2D(x, p)/l);
	else
		w1 = 1.0;
	
	_intersect( y, p, z, x, I );
	l = _length2D(y,I);
	if( l != 0.0 ) 
		w2 = 1.0 - (_length2D(y, p)/l);
	else
		w2 = 1.0;
	
	_intersect( z, p, x, y, I );
	l = _length2D(z,I);
	if( l != 0.0 )
		w3 = 1.0 - (_length2D(z, p)/l);
	else
		w3 = 1.0;
	
	return w1 * f(x, p[0],p[1]) + w2 * f(y, p[0],p[1]) + w3 * f(z, p[0],p[1]);
}
    
    
// Find intersection between AB-> and CD->
void CorrectionSpaceMesh::Triangle::_intersect( osg::Vec3d A, osg::Vec3d B, 
                 osg::Vec3d C, osg::Vec3d D,
                 osg::Vec3d &I )
{
	double m1, m2;
	double d1, d2;
	bool m1i, m2i;  // m1 is not infinite, m2 is not infinte;
	
	m1i = _eqline( A, B, m1, d1 );
	
	m2i = _eqline( C, D, m2, d2 );
	
	// If the first slope (m1) is infinite
	if( m1i == false )
	{
		I[0] = A[0];
		I[1] = m2 * I[0] + d2;
	}
	// if the second slope (m2) is infinite
	else if( m2i == false )
	{
		I[0] = C[0];
		I[1] = m1 * I[0] + d1;
	}
	// Only other case.  In unconstrained space, lines could also be 
	// parallel, but in the context of this definition they cannot
	else
	{
		I[0] = (d2 - d1)/(m1 - m2);
		I[1] = m1 * I[0] + d1;
	}
}

bool CorrectionSpaceMesh::Triangle::_eqline( osg::Vec3d A, osg::Vec3d B, double &m, double &d )
{
	double dx = B[0] - A[0];
	double dy = B[1] - A[1];
	// If slope is infinite (dx == 0), return false, equation is y = A[0];
	if( dx == 0.0 )
		return false;
	
	m = dy/dx;
	
	// d = y  - mx;
	d = A[1] - m * A[0];
	return true;
}

double CorrectionSpaceMesh::Triangle::_length2D( osg::Vec3d A, osg::Vec3d B )
{
	return sqrt( sqr(B[0] - A[0]) + sqr(B[1] - A[1]) );
}
