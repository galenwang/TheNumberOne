#include "common.h"
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgDB/FileNameUtils>
#include <osgDB/WriteFile>
#include <osg/CoordinateSystemNode>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/LineSegment>
#include <osg/NodeVisitor>
#include <osgUtil/IntersectVisitor>
#include <osg/TriangleFunctor>
#include <osgUtil/DelaunayTriangulator>
#include <osgSim/LightPointNode>
#include <osgSim/Sector>
#include <osg/MatrixTransform>

using std::string;

#define	M_PI	3.14159265358979323846			/* 这个值搞错了，导致计算偏差迟迟查不到原因。现予以纠正 2012-01-26 */
#define M_SQRT2	sqrt(2)/2

#define _SQUASH    0.9966471893352525192801545
#define _STRETCH   1.0033640898209764189003079
#define _POLRAD    6356752.3142451794975639668
#define _EQURAD    6378137.0
#define _FLATTENING 298.257223563

#ifndef M_PI
#define SGD_PI 3.14159265358979323846   /* From M_PI under Linux/X86 */
#else
#define SGD_PI M_PI
#endif
#define SGD_180   180.0
#define SGD_DEGREES_TO_RADIANS  (SGD_PI/SGD_180)
#define SGD_RADIANS_TO_DEGREES  (SGD_180/SGD_PI)

#define E2 fabs(1 - _SQUASH*_SQUASH)

#define NEED_ELEVATION  1
#define NO_ELEVATION    0 

#define DEFINE_FOREST_STRING 	"EvergreenForest"
#define DEFINE_FOREST_STRING2	"Forest"
#define DEFINE_LAKE_STRING 		"Lake"
#define DEFINE_OCEAN_STRING 	"Ocean"
#define DEFINE_LIGHTS_STRING 	"Lights"
#define DEFINE_ROAD_STRING 		"Road"
#define DEFINE_RIVER_STRING 	"River"
#define DEFINE_DETAIL_STRING 	"DETAILEDTXT"
#define DEFINE_ISLAND_STRING 	"Island"
#define DEFINE_AIRPORT_STRING 	"Airport"
#define DEFINE_DYNAMICGRASS_STRING 	"DynamicGrass"

/******************************************************************
 *
 *    管理用户SHP文件
 *
 ******************************************************************/
typedef struct _SHPMANAGEMENT_
{
	string fileName;
	int btgID;
	int longi;
	int lati;
	int type;				//0:Point;	1:Line;		2:Polygon;
}ShpManagement;


/******************************************************************
 *
 *    随机数生成基本算法
 *
 ******************************************************************/
#define MT_N 624
#define MT_M 397

typedef struct _MT_
{
	unsigned int array[MT_N]; 
	int index; 
}mt;

void mt_init(mt *mt, unsigned int seed);

/******************************************
 *
 *    地心直角坐标系、经纬高坐标系数据结构
 *
 ******************************************/

typedef struct _Carton_
{
	double x;
	double y;
	double z;
}Carton;

typedef	struct _Geodesy_
{
	double _lon;
	double _lat;
	double _elevation;
}Geodesy;

double deg2rad(double val);

double rad2deg(double val);

double pi();

// given lat1, lon1, lat2, lon2, calculate starting and ending
// az1, az2 and distance (s).  Lat, lon, and azimuth are in degrees.
// distance in meters
static int _geo_inverse_wgs_84( double lat1, double lon1, double lat2,
			double lon2, double *az1, double *az2,
                        double *s );

double distanceM(Geodesy *p1, Geodesy *p2);

void CartToGeod(Carton *cart, Geodesy *geod);

void GeodToCart(Geodesy *geod, Carton *cart);

double distance3Dsquared(Carton *a, Carton *b);

double distance3D(Carton *a, Carton *b);


/******************************************
 *
 *          BUCKET 数据结构
 *
 ******************************************/

/**
 * standard size of a bucket in degrees (1/8 of a degree)
 */
#define SG_BUCKET_SPAN  0.125
#define SG_EPSILON 	0.0000001

/**
 * half of a standard SG_BUCKET_SPAN
 */
#define SG_HALF_BUCKET_SPAN ( 0.5 * SG_BUCKET_SPAN )


// return the horizontal tile span factor based on latitude
static double sg_bucket_span( double l );

typedef struct _Bucket_
{
    short lon;        // longitude index (-180 to 179)
    short lat;        // latitude index (-90 to 89)
    char x;          // x subdivision (0 to 7)
    char y;          // y subdivision (0 to 7)
}Bucket;

long int gen_index(Bucket *B); 

void gen_index_str(Bucket *B, char *idx);

double get_center_lat(Bucket *B); 

// return width of the tile in degrees
double get_width(Bucket *B); 

// return height of the tile in degrees
double get_height(Bucket *B); 

// Build the path name for this bucket
void gen_base_path(Bucket *B, char *path);

double get_center_lon(Bucket *B);

void set_bucket(Bucket *B, double dlon, double dlat );

void SGBucket(Bucket *B, const long int bindex);

/**********************************************************
 *
 *  
 *        BTG 数据结构
 *
 *
 **********************************************************/
typedef struct _TERPOINT_
{
	double px;
	double py;				
	double pz;				//实际坐标   需要计算转换过来的
	double lon;
	double lat;				
	double elev;			//经纬度坐标，需要计算得来
	short  nx;
	short  ny;
	short  nz;				//法向量， 
}TerPoint;

typedef struct _TERTEXTCOORD_
{
	float  fU;
	float  fV;
}TerTextcoord;

typedef struct _TRIVERTEX_
{
	short P;
	short T;
	short N;
}TriVertex;

typedef struct _TRIELEMENT_
{
	int iSum;
	int isBound;
	std::vector<TriVertex> mVertex;
	short minP;
}TriElement;

typedef struct _TRITYPE_
{
	int iElems;
	char cMtype[32];
	int cObjType;
	int cPropVal;
	std::vector<TriElement> mElement;
}TriType;

int ReadBTGToDataStruct(char *fname);

int WriteDataStructToBTG(char *fname);

int WriteDataStructToBTG_River(char *fname);

int WriteDataStructToBTG_OnlyCulture(char *fname);

int WriteDataStructToBTG_CulturePlusIsland(char *fname);

/******************************************
 *
 *          SimEdges 
 *
 ******************************************/
struct SimEdge
{
	osg::Vec3d	A, B;
	int			x, y;
	int			z;
	bool		keep;

    bool operator == (SimEdge &rhs )
    {
        if( (A == rhs.A && B == rhs.B) || (A == rhs.B && B == rhs.A))
            return true;
        else
            return false;
    }

    bool has( osg::Vec3d v )
    {
        if( v == A || v == B ) return true;
        return false;
    }
    
    int  get( osg::Vec3d v )
    {
    	if( v == A ) return z;
    	else		 return -1;
    }
};

struct SimMakeEdgeList
{
    public:
        SimMakeEdgeList():edges(0L) {}

        void setEdgeVector( std::vector<SimEdge> *e )
        {
            edges = e;
        }

        void operator() (const osg::Vec3d& v1,const osg::Vec3d& v2,const osg::Vec3d& v3, bool flag) const
        {
            if( edges != 0L )
            {
				SimEdge se; 
				se.A = v1; se.B = v2; se.keep = true; edges->push_back(se);
				se.A = v2; se.B = v3; se.keep = true; edges->push_back(se);
				se.A = v3; se.B = v1; se.keep = true; edges->push_back(se);
            }
        }

    private:
        std::vector<SimEdge> * edges;
};

struct GeoEdge
{
	Geodesy		A, 	B;
	osg::Vec3d  Ac, Bc;
	double		dist;
	
    void distance()
    {
    	dist = distanceM(&A, &B);
    }    
};

/******************************************
 *
 *          IndEdges 
 *
 ******************************************/
struct IndEdge
{
	int  A, B;
	bool keep;

    bool operator == (IndEdge &rhs )
    {
        if( (A == rhs.A && B == rhs.B) || (A == rhs.B && B == rhs.A))
            return true;
        else
            return false;
    }

    bool has( int v )
    {
        if( v == A || v == B ) return true;
        return false;
    }
};

struct IndMakeEdgeList
{
    public:
        IndMakeEdgeList():edges(0L) {}

        void setEdgeVector( std::vector<IndEdge> *e )
        {
            edges = e;
        }

        void operator() (const int v1,const int v2,const int v3) const
        {
            if( edges != 0L )
            {
				IndEdge se; 
				se.A = v1; se.B = v2; se.keep = true; edges->push_back(se);
				se.A = v2; se.B = v3; se.keep = true; edges->push_back(se);
				se.A = v3; se.B = v1; se.keep = true; edges->push_back(se);
            }
        }

    private:
        std::vector<IndEdge> * edges;
};


/******************************************
 *
 *          FoundNode 
 *
 ******************************************/

typedef struct _FOUNDNODE_
{
	string 	name;
	int    	gmIdx;
	double  maxX;
	double  maxY;
	double  maxZ;
	double  minX;
	double  minY;
	double  minZ;
}FoundNode;

typedef struct _MAXMINLONLAT_
{
	string mShpName;
	double maxLong;
	double maxLati;
	double minLong;
	double minLati;
}MaxMinLonLat;

void InitMMLLEx(MaxMinLonLat *m, string fn);
void SetMMLLEx(MaxMinLonLat *m, double lon, double lat);
/******************************************
 *
 *          Point的文化信息 
 *
 ******************************************/
typedef struct _SIMPOINT_
{
	string name;
	double lon;
	double lat;
	double elev;
	double x;
	double y;
	double z;
	string tile;
	int index;
}SimPoint;

typedef struct _LIGHTVERTEX_
{
	osg::Vec3d	position;
	osg::Vec3d	direction;
	osg::Vec4	color;	
	int 		index;
}LightVertex;

typedef struct _LIGHTNODE_
{
	string 		name;
	std::vector<LightVertex> vlist;
}LightNode;

typedef struct _COLORNAME_
{
	string      name;
	osg::Vec4	color;	
}ColorName;

/******************************************
 *
 *          Function 
 *
 ******************************************/


int ChangeVertexInDataStru(Geodesy *source);

int FindSkirtVertex(float elev);

/*从DataStru去掉裙边的面*/
int RemoveSkirtsFromDataStru();

int CalculatePerimeter( osg::Vec3Array &coords );

osg::ref_ptr<osg::Geode> GetAllLineLoopFromSkirtEdge();

osg::ref_ptr<osg::Group> CreateGroupFromDataStru();

int VerifyTheModels(char *fn, char *out);

int ReadFltAndVistIt(char *fname, osg::ref_ptr<osg::Geode> gnode);

void ReadIVEAndChangeVertexToGeodesy(char *fname);
void ReadIVEAndComputePerimeter(char *fname);
void ReadIVEAndDeleteLine(char *fname);
void ReadIVEAndConnectLineLoop(char *fname);                        
