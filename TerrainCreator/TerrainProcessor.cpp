#include "common.h"
#include <osg/CoordinateSystemNode>
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osg/MatrixTransform>
#include <osg/NodeVisitor>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/PagedLOD>
#include <osg/Switch>
#include <osg/Texture>
#include <osg/Texture1D>
#include <osg/Texture2D>
#include <osg/Texture3D>
#include <osgFX/MultiTextureControl>
#include <osgUtil/DelaunayTriangulator>
#include <osgUtil/SmoothingVisitor>
#include <osgUtil/IntersectVisitor>
#include <osgUtil/Tessellator>
#include <osg/AlphaFunc>
#include <iostream>
#include <assert.h>
#include <math.h>
#include "TerrainProcessingUtils.h"
#include "TexCoordPreserver.h"
#include "CorrectionSpaceMesh.h"
#include "bucket.h"


#define ADD_SKIRT
//#define MAKE_HIGH_SKIRT
extern double getArea(osg::Vec3 v0, osg::Vec3 v1, osg::Vec3 v2);
extern osg::Vec3 cross(const osg::Vec3 v1, const osg::Vec3 v2);
extern double distanceM(Geodesy *p1, Geodesy *p2);
extern int ReadBTGToDataStruct(char *fname);             
extern osg::ref_ptr<osg::Group> createGroupFromDataStru();
extern std::vector< osg::ref_ptr<CorrectionSpaceMesh> > csmList;
extern char tempDirectory[_MAX_PATH];
extern char debugFile[_MAX_PATH];
/* 把一个字符串按照下划线或者中线分成各个部分。*/
extern std::vector<std::string> materialParts;
extern int splitMaterialName(std::string matname);
int CurrentTerrainIsCut = 0;
int IsTopTerrain;
int CurrentTerrainLodLevel = 0;
int UpdateTheElevationInTerrain = 0;	//如果更新了地景块内某个点的高程。
osg::ref_ptr<osg::Group> allCultures;
double elevationOfLake = -9999.99;		//记录当前的湖面海拔高度。
int totalRoad = 0;
int totalLake = 0;
int totalIsland = 0;
#define SKIRT_HEIGHT 200.0
#define UNITDISTANCE 10.0
#define AREA_LIMIT 1000.0
double m_elevOfAirport;

std::string cultures[8] = 
{ 
	"Airport",
	"river",
	"River",
	"Lake",
	"Stream",
	"Road",
	"Island",
	"DynamicGrass"
};
int culture_cc = 8;

typedef struct _INTERSDIFF_
{
	osg::Vec3 intrs;
	int diff;
	int border;
}IntersDiff;

std::vector<std::string> detailedtxt;
int detailedtxt_cc = 0;

/*  解析线的名称，以下划线或者中线作为分段标志，将名称分为几段，*/
std::vector<std::string> NameParts;
void SplitName(std::string matname)
{
	char matpart[128], mat[256];
	int j = 0;
	
	/* 把matname以下划线作为分隔符，分成若干小段，保存在nameParts里面 */
	sprintf(mat, "%s", matname.c_str());
	NameParts.clear();
	for(unsigned int i = 0; i<strlen(mat); i++)
	{
		if((mat[i] != '_')&&(mat[i] != '-')) { matpart[j++] = mat[i];	}
		if((mat[i] == '_')||(mat[i] == '-')) { matpart[j] = 0;	NameParts.push_back(matpart);	j = 0; }
	}
	matpart[j] = 0;	NameParts.push_back(matpart);	j = 0;	
}

/* 获取三个点的法线。*/
osg::Vec3 getNormal(osg::Vec3 v0, osg::Vec3 v1, osg::Vec3 v2)
{
	osg::Vec3 normal = cross(v1 - v0, v2 - v0);
	return normal;
}

/* 得到一个线段长度值 */
double distanceV(osg::Vec3 v0, osg::Vec3 v1)
{
	return sqrt( ((v0.x() - v1.x()) * (v0.x() - v1.x())) +  ((v0.y() - v1.y()) * (v0.y() - v1.y())) +  ((v0.z() - v1.z()) * (v0.z() - v1.z())) );
}

/* 检查一个字符串是否是细节纹理字符串。 */
bool isDetailedtxt(std::string name)
{
	bool isDetailedtxt = false;
	for (int i=0;i<detailedtxt_cc;i++)
	{
		if (name.substr(0, strlen(detailedtxt[i].c_str())) ==detailedtxt[i])
		{
			isDetailedtxt = true;
			break;
		}
	}

	return isDetailedtxt;
}

std::vector<osg::ref_ptr<osg::Geode>> detailed_gdlist;
std::vector<osg::ref_ptr<osg::Geode>> dynamicGrass_gdlist;

/* 从一个Node里面找到最顶层的特定类型 */
template<class T>
class FindTopMostNodeOfTypeVisitor : public osg::NodeVisitor
{
public:
    FindTopMostNodeOfTypeVisitor():
        osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN),
        _foundNode(0)
    {}

    void apply(osg::Node& node)
    {
        T* result = dynamic_cast<T*>(&node);
        if (result)
        {
            _foundNode = result;
        }
        else
        {
            traverse(node);
        }
    }

    T* _foundNode;
};

template<class T>
T* findTopMostNodeOfType(osg::Node* node)
{
    if (!node) return 0;

    FindTopMostNodeOfTypeVisitor<T> fnotv;
    node->accept(fnotv);

    return fnotv._foundNode;
}

class GeodeInfo
{
public:
    GeodeInfo( float d, osg::Geode* g )
        : _distance( d ), _geode( g ) {}

    float _distance;
    osg::ref_ptr<osg::Geode> _geode;
};
typedef std::vector<GeodeInfo> GeodeList;

/* 遍历Node，查找文化信息或者细节纹理的node部分。保存在_features里面，其余的保存在_nonFeatures里面 */
class FindFeatures: public osg::NodeVisitor
{
    public:

        FindFeatures()
          : osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN ),
            _switchInDistance( 1000.f )
        {}

        virtual void apply( osg::Geode &geode )
        {
            GeodeInfo gi( _switchInDistance, &geode );

            // For Now look at the name.
            std::string name = geode.getName();
            if( _isFeature( geode ) )
                _features.push_back( gi );
            else 
                _nonFeatures.push_back( gi );

            traverse( geode );
        }

        virtual void apply( osg::LOD& lod )
        {
            float lastDistance = _switchInDistance;

            int i;
            for (i=0; i<(int)lod.getNumChildren(); i++)
            {
                if( lod.getRangeMode() == osg::LOD::DISTANCE_FROM_EYE_POINT )
                {
                    _switchInDistance = lod.getMinRange( i );
                    if (_switchInDistance < 1000.f) _switchInDistance = 1000.f;
                }
                else
                {
                    // Range mode is pixel size. Hm. Not sure what to do here.
                    _switchInDistance = 1000.f;
                    osg::notify( osg::WARN ) << 
                        "Terrain Deformation Software Loader: Terrain Processor:\n" <<
                        "  Scene graph contains LOD Node with unsupported range mode." << std::endl;
                }
                lod.getChild( i )->accept( *this );
            }

            _switchInDistance = lastDistance;
        }

        const GeodeList &getFeatures() const { return _features; }
        const GeodeList &getNonFeatures() const { return _nonFeatures; }

    private:
        bool _isFeature( const osg::Node& node )
        {
					// Look for feature as node name
					
					std::string name = node.getName();
					for (int i=0; i<culture_cc; i++)
					{
						if (name.substr(0, strlen(cultures[i].c_str())) == cultures[i] )
							return true;
					}
					for (int i=0; i<detailedtxt_cc; i++)
					{
						if (name.substr(0, strlen(detailedtxt[i].c_str())) == detailedtxt[i])
							return true;
					}
					
					// Look for feature in description string. (OpenFlight stores
					//   comment record text in Node's DescriptionList.)
					int idx = node.getNumDescriptions();
					while (idx--)
					{
						const std::string& desc = node.getDescription( idx );
						
						for (int i=0; i<culture_cc; i++)
						{
							if (desc.find( cultures[i] ) != desc.npos)
								return true;
						}
						for (int i=0; i<detailedtxt_cc; i++)
						{
							if (desc.find( detailedtxt[i] ) != desc.npos)
								return true;
						}
					}
					
					return false;
        }

        GeodeList _features;
        GeodeList _nonFeatures;

        float _switchInDistance;
};

/* 德罗尼约束 */
class ArealConstraint: public  osgUtil::DelaunayConstraint 
{ 
public:
	int getinteriorTrisSize() { return _interiorTris.size(); }

	void handleOverlaps1(void)
	{
		// use tessellator to interpolate crossing vertices.
		osg::ref_ptr<osgUtil::Tessellator> tscx=new osgUtil::Tessellator; // this assembles all the constraints
		tscx->setTessellationType(osgUtil::Tessellator::TESS_TYPE_GEOMETRY);
		tscx->setBoundaryOnly(true);
		tscx->setWindingType( osgUtil::Tessellator::TESS_WINDING_ABS_GEQ_TWO); 
	    
		tscx->retessellatePolygons(*this); // find all edges
	}
};

std::vector<FindValidVertex> terrainVertices;
std::vector<FindValidVertex> terrainBorderVertices;
std::vector<FindValidVertex> terrainInnerVertices;
std::vector<FindValidVertex> allBorderVertices;
std::vector<FindValidVertex> faceVertice;
std::vector<FindValidVertex> skirtVertics;
std::vector<FindValidVertex> validVertices;
std::vector<FindValidVertex> validBorderVertices;
std::vector<FindValidVertex> intersectlist;
std::vector<FindBorderVertex> checkBorderVertices;
std::vector<osg::Vec3> invalidVertices;
std::vector<osg::Vec3> innerVerticesInterWithCulture;
std::vector<unsigned int> faceElement;
std::vector<unsigned int> skirtElement;
std::vector<unsigned int> validElement;
std::vector<osg::ref_ptr<osg::Vec3Array>> perimeterList;
osg::ref_ptr<osg::Vec3Array> terrainBorderLineLoop;
osg::ref_ptr<osg::Vec3Array> terrainBoundLoop;
std::vector<osg::ref_ptr<ArealConstraint> > dclist;
std::vector<osg::ref_ptr<ArealConstraint> > dcDetail;
std::vector<osg::ref_ptr<ArealConstraint> > dcLakes;
std::vector<osg::ref_ptr<osg::Vec3Array> > multiDc;
osg::Geometry *pCurTerrain;
osg::ref_ptr<osg::Geode> pCurrTerrainGeode;
bool AllBorderVertexIntersection;				//所有边界的点都跟文化信息相交？
bool AllInnerVertexIntersection;				//所有内部的点都跟文化信息相交？
bool AllInnerVertexNoIntersection;				//所有内部的点都跟文化信息不相交？
int IntersectCt;												//记录地景面和文化信息有多少个交点的全局变量。
int subIntersectCt;											//
bool needDoSkirt;												//决定一条外围环是否可以制作裙边。
double skirtLength;											//裙边的高度 m
osg::BoundingSphere terrainBS;					//地景块的包围球
std::vector< osg::ref_ptr<CorrectionSpaceMesh> > validCSM;

std::vector<ProjectionMap> totalVertices;
std::vector<ProjectionMap> borderPMs;
std::vector<ProjectionMap> featurePMs;
std::vector<ProjectionMap> IntersectVertices;		//记录所有的交点，如果下一个交点的投影点，同前面某个交点投影点相同，则直接使用那个交点的空间点信息，避免重新计算产生错层。
std::vector<ProjectionMap> roadTexCoordSave;		//记录所有路顶点的纹理坐标

#ifdef ADD_SKIRT
void CreateSkirt(osg::ref_ptr<osg::Vec3Array> validCoords, osg::ref_ptr<osg::DrawElementsUInt> dui, int mode);
#endif

/* 线段数据结构 */
typedef struct _LINESEG_
{
	osg::Vec3 pS;
	osg::Vec3 pE;
	int idxS;
	int idxE;
	int validS;
	int validE;
	int border;
}LineSeg;

/* Add by Galen 2013-02-21 */
/* 下面数据结构保存顶点到湖边段的距离。*/
typedef struct _VERTEXDISTANCE_
{
	osg::Vec3 pVF;
	double distance;
}VertexDistance;
std::vector<VertexDistance> vertexNearLake;
#define MINDISTANCETOEDGE 20.0
std::vector<TerrainProcessingUtils::Edge> allLakeEdges;
/* Add by Galen 2013-02-21 End */

/* Add by Galen 2013-02-22 */
/* 下面数据结构保存已经边段化的湖的文化信息 */
typedef struct _LAKEEDGEITEM_
{
	int numOfVertex;
	osg::Vec3 firstVertex;
}LakeEdgeItem;
std::vector<LakeEdgeItem> allEdgedLakes;
/* Add by Galen 2013-02-22 End */

FindValidVertex conner[4];
FindValidVertex connerEx[4];						//根据conner四个点求得，但比这四个点略微扩大一圈的点。
LineSeg border[4];
std::vector<LineSeg> borderIntersectLS;				//记录所有的边段，原始边段，包括与文化信息相交(两个端点都相交)、半相交(只有一个端点相交)和不相交(两个端点都不相交)的所有边段；
std::vector<LineSeg> borderHalfIntersectLS;		//记录所有的半相交(只有一个端点相交)边段；
std::vector<LineSeg> validBorderSegment;			//记录所有有效的边段，包括计算了插入点之后的有效边；

int f_idx;

osg::ref_ptr<osg::Vec3Array> FindValidVertexInDC(osg::ref_ptr<ArealConstraint> dc);
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToTranglesFace(osg::ref_ptr<osg::Vec3Array> gvx);

/* 检查一个点是否在FindValidVertex点集合内 */
int checkElementUI(std::vector<FindValidVertex> fvvlist, int ia, int ib, int ic)
{
	int iFlag = 1;
	for( std::vector<FindValidVertex>::iterator f= fvvlist.begin(); f != fvvlist.end(); f++ )
	{
		if (ia == f->idx)  { iFlag = 0; break; }
		if (ib == f->idx)  { iFlag = 0; break; }
		if (ic == f->idx)  { iFlag = 0; break; }
	}
	return iFlag;
}



/* 根据conner来计算connerEx。*/
void calcConnerEx()
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 v, r;
	double maxLon, maxLat, minLon, minLat;
	maxLon = maxLat = -999999.999;	minLon = minLat = 999999.999;
	ProjectionMap cnEx[4], *pm = &cnEx[0];
	
	/* 先求出四个角的经纬度，及相关极值。*/
	for(unsigned int i=0; i<4; i++)
	{
		v = conner[i].pV;
		cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
		CartToGeod(&cat, &geo);	pm->lon = geo._lon; pm->lat = geo._lat; pm->elev = geo._elevation;

		if(geo._lon > maxLon) maxLon = geo._lon;
		if(geo._lon < minLon) minLon = geo._lon;
		if(geo._lat > maxLat) maxLat = geo._lat;
		if(geo._lat < minLat) minLat = geo._lat;
		pm++;
		
		connerEx[i].pT = conner[i].pT;	connerEx[i].idx = conner[i].idx;	connerEx[i].newIdx = conner[i].newIdx; connerEx[i].isValid = conner[i].isValid;
	}
	
	/* 然后写ConnerEx */
	pm = &cnEx[0];
	for(unsigned int i=0; i<4; i++)
	{
		if(fabs(pm->lon - maxLon) < 0.00001)  pm->lon += 0.001;
		if(fabs(pm->lon - minLon) < 0.00001)  pm->lon -= 0.001;
		if(fabs(pm->lat - maxLat) < 0.00001)  pm->lat += 0.001;
		if(fabs(pm->lat - minLat) < 0.00001)  pm->lat -= 0.001;
		geo._lon = pm->lon; geo._lat = pm->lat; geo._elevation = pm->elev;	GeodToCart(&geo, &cat);	pm->pV.set(cat.x, cat.y, cat.z);
		geo._lon = pm->lon; geo._lat = pm->lat; geo._elevation = 0.0;		GeodToCart(&geo, &cat);	pm->pVF.set(cat.x, cat.y, cat.z);
		
		connerEx[i].pV = pm->pV;	connerEx[i].pVF = pm->pVF;
		pm++;
	}
	
	/* 把这几个点加入到无效点集合里面 */
	invalidVertices.push_back(connerEx[0].pVF);
	invalidVertices.push_back(connerEx[1].pVF);
	invalidVertices.push_back(connerEx[2].pVF);
	invalidVertices.push_back(connerEx[3].pVF);
}

/* 输入一个投影点，找到该点对应的3D空间点。*/
osg::Vec3 getSpaceVertex(osg::Vec3 v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 r = v;
	unsigned int i, iGet = 0;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	
	for(i=0; i<totalVertices.size(); i++)
	{
		if(	(fabs(geo._lon - totalVertices[i].lon) < 0.00001) &&
			(fabs(geo._lat - totalVertices[i].lat) < 0.00001) )
		{	iGet = 1;	break;	}
	}
	if(iGet)	r = totalVertices[i].pV;
	return r;
}

/* 输入一个投影点，找到该点对应的3D空间点。*/
osg::Vec3 getSpaceVertexEx(osg::Vec3 v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 r = v;
	unsigned int i, iGet = 0;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	
	for(i=0; i<totalVertices.size(); i++)
	{
		if(	(fabs(geo._lon - totalVertices[i].lon) < 0.00001) &&
			(fabs(geo._lat - totalVertices[i].lat) < 0.00001) )
		{	iGet = 1;	break;	}
	}
	
	if(iGet)
	{
		/* 高度值相加 */
		geo._elevation +=  totalVertices[i].elev;
		GeodToCart(&geo, &cat);		
		r.x() = cat.x;	r.y() = cat.y;	r.z() = cat.z;
	}
	return r;
}


/* 输入一个Vertex， 输出转换成经纬高坐标系后，将高度值改为0.0的点，即求出该点在海平面的映射点。*/
osg::Vec3 getProjectionVertex(osg::Vec3 v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 r;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	geo._elevation = 0.0; GeodToCart(&geo, &cat);
	r.x() = cat.x; r.y() = cat.y; r.z() = cat.z;
	
	return r;
}

/* 输入一个Vertex， 计算出ProjectionMap结构*/
void getProjectionMapByVertex(ProjectionMap *pm, osg::Vec3 v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 r;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	pm->elev = geo._elevation;	geo._elevation = 0.0; GeodToCart(&geo, &cat);
	r.x() = cat.x; r.y() = cat.y; r.z() = cat.z;
	
	pm->pV = v;	pm->pVF = r; pm->lon = geo._lon; pm->lat = geo._lat; 
}

/* 根据当前块是否有湖，来调整当前块各个点的高度。原理和调整csm面高度相同。*/
void TuneElevationOfLakeSurface(osg::ref_ptr<osg::Geode> lakeSurface, osg::ref_ptr<osg::Vec3Array> coords)
{
	if(lakeSurface->getNumDrawables() > 0)
	{
		Geodesy geo;
		Carton  cat;
		osg::Vec3 v;
		
		/* Add by Galen 2013-02-21 处理湖边的点。湖边点要尽量平滑过度一下，避免太突兀。*/
		if(vertexNearLake.size() > 0)
		{
			for(std::vector<VertexDistance>::iterator i = vertexNearLake.begin(); i != vertexNearLake.end(); ++i)
			{
				/* 首先要根据这个点，从faceVertice里面找到原始高度 */
				double currVerHeight;
				for( std::vector<FindValidVertex>::iterator f = faceVertice.begin(); f != faceVertice.end(); f++ )
				{
					if(f->pVF == i->pVF)
					{
						cat.x = f->pV.x(); cat.y = f->pV.y(); cat.z = f->pV.z();
						CartToGeod(&cat, &geo);	currVerHeight = geo._elevation;
						
						/* 根据线性法则计算调整后的湖边点高度，并写回 faceVertice 里面 */
						geo._elevation = elevationOfLake + (currVerHeight - elevationOfLake) * ( i->distance / MINDISTANCETOEDGE);
						GeodToCart(&geo, &cat);
						v.x() = cat.x; v.y() = cat.y; v.z() = cat.z;
						f->pV = v;
						terrainVertices[f->idx].pV = f->pV;						 
						
						break;
					}
				} // for f
			}  // for i
		}
		/* Add by Galen 2013-02-21 End */

		/* 处理地景面上的点 */
		for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		{
			osg::Vec3 isect;
			v = f->pVF;
			if( TerrainProcessingUtils::getIntersect( lakeSurface, v, isect ) != false )
			{
				cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
				CartToGeod(&cat, &geo);	geo._elevation = elevationOfLake; GeodToCart(&geo, &cat);
				v.x() = cat.x; v.y() = cat.y; v.z() = cat.z;
				f->pV = v;
				terrainVertices[f->idx].pV = f->pV;
			}
		}
		
		/* 处理裙边上的点 */
		for( std::vector<FindValidVertex>::iterator f= skirtVertics.begin(); f != skirtVertics.end(); f++ )
		{
			osg::Vec3 isect;
			v = f->pVF;
			if( TerrainProcessingUtils::getIntersect( lakeSurface, v, isect ) != false )
			{
				cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
				CartToGeod(&cat, &geo);	geo._elevation = elevationOfLake - skirtLength; GeodToCart(&geo, &cat);
				v.x() = cat.x; v.y() = cat.y; v.z() = cat.z;
				f->pV = v;
				terrainVertices[f->idx].pV = f->pV;
			}
		}
		
		/* 把有高度变化的点写回到Geometry里面。*/
		int index = 0;
		for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
		{
			*r = (terrainVertices[index++].pV);
		}		
	}
}


/* 已知一组点，类型为FindValidVertex，计算pV对应的pVF，即这个点转换为经纬高之后，对应高度为0的那个点。*/
void getFlatVertex(std::vector<FindValidVertex> *plist)
{
	for( std::vector<FindValidVertex>::iterator p= plist->begin(); p != plist->end(); p++ )
	{
		osg::Vec3 pF = getProjectionVertex(p->pV);
		p->pVF = pF; 
	}
}

/* 已知一组点，类型为FindValidVertex，把里面所有投影点pVF，组成一个osg::Vec3Array */
osg::ref_ptr<osg::Vec3Array> getCoordsFlat(std::vector<FindValidVertex> plist)
{
	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator p= plist.begin(); p != plist.end(); p++ )
	{
		vx->push_back(p->pVF);
	}
	return vx.release();	
}

/* 把当前feature的各个Geometry里面所有的Vertex都改成经纬高坐标里面高度为0的对应点，并存储到featurePMs里面。 */
void makeGeometryFlat(osg::ref_ptr<osg::Geometry> fg)
{
	osg::ref_ptr<osg::Vec3Array> vx = dynamic_cast<osg::Vec3Array *>(fg->getVertexArray());
	for( osg::Vec3Array::iterator r = vx->begin(); r != vx->end(); r++ )
	{
		/* 增加文化信息点的映射, 便于后面查找. */
		ProjectionMap pm;
		getProjectionMapByVertex(&pm, *r);
		featurePMs.push_back(pm);	

		osg::Vec3 pF = pm.pVF;
		*r = pF;
	}
}

/* 根据一个feature点，从featurePMs里面找到对应项，并把这个对应项保存到 totalVertices 里面。*/
osg::Vec3 findPMByFeatureVertex(osg::Vec3 vx, int rettype, string pName)
{
	for( std::vector<ProjectionMap>::iterator p = featurePMs.begin(); p != featurePMs.end(); p++ )
	{
		if((vx == p->pV) || (vx == p->pVF))	
		{
			/* 只有当文化信息名字为路、岛或者细节纹理面的时候，才调整文化信息点的实际高度同地形面高度相同。*/
			/* 如果是其他情形，如湖、机场外围环，仍按原来的文化信息点高度计算。*/
			if(		(pName.substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING) || 
					(pName.substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) || 
					(isDetailedtxt(pName)) )
			{
				/* 这里增加一段代码，重新计算文化信息点的高度 2012-11-25 */
				/* 计算方法是同当前这个地形面做碰撞检测，修改featurePMs里面的原始值。*/
				/* 创建相交线段，是一根在经度lon和纬度lat的垂直的线段。用这根线段与地景做相交。*/
				Geodesy stg; stg._lon = p->lon; stg._lat = p->lat; stg._elevation =  9999999.0;
				Carton  stc;
				GeodToCart(&stg, &stc);
				osg::Vec3d start; start.set(stc.x, stc.y, stc.z); 
				Geodesy eng; eng._lon = p->lon; eng._lat = p->lat; eng._elevation =  -9999999.0;
				Carton  enc;
				GeodToCart(&eng, &enc);
				osg::Vec3d end;   end.set(enc.x, enc.y, enc.z); 
	
				osgUtil::IntersectVisitor intersectVisitor; 
				osg::LineSegment* lineSegment = new osg::LineSegment(start, end); 
				intersectVisitor.addLineSegment(lineSegment); 
				pCurrTerrainGeode->accept(intersectVisitor); 
				
				/* 获取相交结果 Get intersections */
				bool hits = intersectVisitor.hits(); 
				if (hits) 
				{ 
					int nHits = intersectVisitor.getNumHits(lineSegment); 
					double alt = -999.0; 
					osg::Vec3 curV;
					for (int i = 0; i < nHits; ++i) 
					{ 
						const osgUtil::Hit& hit = intersectVisitor.getHitList(lineSegment)[i]; 
						osg::Vec3d point; 
						point = hit.getWorldIntersectPoint(); 
						Carton pntc; pntc.x = point.x(); pntc.y = point.y(); pntc.z = point.z();
						Geodesy pntg;
						CartToGeod(&pntc, &pntg);
						double elevation = pntg._elevation; 
						if (alt < elevation) 
						{ 
							alt = elevation; 
							curV.set(pntc.x, pntc.y, pntc.z);
						} 
					}
					
					p->pV = curV;
					p->elev = alt;
				}
				/* 增加完毕。*/
			}else if (pName.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			/* 当文化信息名字为湖的时候， 也按原来的文化信息点高度计算。 */ 
			{
			}

			totalVertices.push_back(*p);
			if(rettype == 0) 	return p->pVF; else return p->pV;
		}
	}
	return vx;
}

/* 把一个交点 插入到 IntersectVertices 里面 */
void addVertexIntoIntersect(osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, v);
	IntersectVertices.push_back(fvv);	
}

/* 从 IntersectVertices 找一个相同的点，找不到就返回-1, 找到就返回下标 */
int FindVertexFromIntersect(osg::Vec3 vF)
{
	if(IntersectVertices.size() == 0) return -1;
	for(unsigned int i=0; i<IntersectVertices.size(); i++)
	{
		float dx = fabs(IntersectVertices[i].pVF.x() - vF.x());
		float dy = fabs(IntersectVertices[i].pVF.y() - vF.y());
		float dz = fabs(IntersectVertices[i].pVF.z() - vF.z());
		if((dx <= 1.0) && (dy <= 1.0) && (dz <= 1.0))
		{
			return i;
		}
	}
	return -1;
}



/* 把一个顶点 插入到 totalVertices 里面 */
void addVertexIntoTotal(osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, v);
	totalVertices.push_back(fvv);	
}

/* 把一个顶点 插入到 totalVertices 里面, 先查找这个顶点在 totalVertices 里面是否有数据 */
void addVertexIntoTotalFirstSearth(osg::Vec3 v)
{
	Geodesy geo, geo2;
	Carton  cat;
	osg::Vec3 r = v;
	unsigned int minIdx, i, iGet = 0;	//minIdx 代表距离这个点最近的那个已知高度的顶点。
	double minSum = 99999.99;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	
	for(i=0; i<totalVertices.size(); i++)
	{
		double deltaLon = fabs(geo._lon - totalVertices[i].lon);
		double deltaLat = fabs(geo._lat - totalVertices[i].lat);
		if(	(deltaLon < 0.00001) &&	(deltaLat < 0.00001) )
		{	iGet = 1;	break;	}

		double deltaSum = deltaLon + deltaLat;
		if(deltaSum < minSum)	{	minSum = deltaSum; minIdx = i;	}
	}
	if(iGet == 0)
	{
		ProjectionMap fvv;
		getProjectionMapByVertex(&fvv, v);
		fvv.pVF = v;
		
		osg::Vec3 vx = totalVertices[minIdx].pV;
		cat.x = vx.x(); cat.y = vx.y(); cat.z = vx.z();
		CartToGeod(&cat, &geo2);	geo._elevation = geo2._elevation;
		GeodToCart(&geo, &cat);		vx.x() = cat.x;	vx.y() = cat.y;	vx.z() = cat.z;
		fvv.pV = vx;
		totalVertices.push_back(fvv);
	}
}

/* Add by Galen 2013-02-18 */
/* 把一个顶点 插入到 totalVertices 里面, 这个顶点的空间高度值为给定的数据 */
void addVertexIntoTotalWithElevation(osg::Vec3 v, float elev)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 r = v;
	unsigned int i, iGet = 0;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	
	for(i=0; i<totalVertices.size(); i++)
	{
		double deltaLon = fabs(geo._lon - totalVertices[i].lon);
		double deltaLat = fabs(geo._lat - totalVertices[i].lat);
		if(	(deltaLon < 0.00001) &&	(deltaLat < 0.00001) )
		{	iGet = 1;	break;	}
	}
	if(iGet == 0)
	{
		ProjectionMap fvv;
		getProjectionMapByVertex(&fvv, v);
		fvv.pVF = v;
		
		osg::Vec3 vx;
		geo._elevation = elev;
		GeodToCart(&geo, &cat);		vx.x() = cat.x;	vx.y() = cat.y;	vx.z() = cat.z;
		fvv.pV = vx;	fvv.elev = elev;
		totalVertices.push_back(fvv);
	}
}
/* Add End */


/* 把一个顶点对插入到 totalVertices 里面 */
void addVertexPairIntoTotal(osg::Vec3 vF, osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, vF);
	fvv.pV = v; fvv.pVF = vF;
	totalVertices.push_back(fvv);	
}


/* 把所有的顶点信息，包括地形点和文化信息点，都集中到一个集合里面，并且计算出每个点的对应海平面点，结果存放在totalVertices里面 */
void getTotalFaceVertex()
{
	/* 先处理地景面上的点 */
	for( std::vector<FindValidVertex>::iterator p= faceVertice.begin(); p != faceVertice.end(); p++ )
	{
		addVertexIntoTotal(p->pV);
	}
}


/* 给定空间两个点，求出它们映射到二维平面上的斜率和截度  */
bool eqline( osg::Vec3d A, osg::Vec3d B, double &m, double &d )
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

/* 将直角坐标系转换为经纬度坐标系 */
osg::Vec3d cTog(osg::Vec3 v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 Vo;

	cat.x = v.x(); cat.y = v.y(); cat.z = v.z();
	CartToGeod(&cat, &geo);	
	Vo.x() = geo._lon;	Vo.y() = geo._lat;	Vo.z() = geo._elevation;
	
	return Vo;
}

/* 将经纬度坐标系转换为直角坐标系 */
osg::Vec3 gToc(osg::Vec3d v)
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 Vo;

	geo._lon = v.x(); geo._lat = v.y(); geo._elevation = v.z();
	GeodToCart(&geo, &cat);	
	Vo.x() = cat.x;	Vo.y() = cat.y;	Vo.z() = cat.z;
	
	return Vo;
}


/* 求出空间两条直线 AB 和 CD 之间的相交情况 */
// Find intersection between AB-> and CD->
void intersect(  osg::Vec3 A, osg::Vec3 B, 
                 osg::Vec3 C, osg::Vec3 D,
                 osg::Vec3 &I )
{
    double m1, m2;
    double d1, d2;
    bool m1i, m2i;  // m1 is not infinite, m2 is not infinite;

#if 1
	/* 先把直角坐标系转换为经纬度坐标系 */
	osg::Vec3d gA, gB, gC, gD, gI;
	gA = cTog(A);	gB = cTog(B);	gC = cTog(C);	gD = cTog(D);

	/* 计算出经纬度格式的交点 */
    m1i = eqline( gA, gB, m1, d1 );
    m2i = eqline( gC, gD, m2, d2 );

    // If the first slope (m1) is infinite
	if( m1i == false )
    {
        gI[0] = gA[0];
        gI[1] = m2 * gI[0] + d2;
    }
    // if the second slope (m2) is infinite
    else if( m2i == false )
    {
        gI[0] = gC[0];
        gI[1] = m1 * gI[0] + d1;
    }
    // Only other case.  In unconstrained space, lines could also be 
    // parallel, but in the context of this definition they cannot
    else
    {
        gI[0] = (d2 - d1)/(m1 - m2);
        gI[1] = m1 * gI[0] + d1;
    }
    
    gI[2] = 0.0;
    I = gToc(gI);
	return;
#else
    m1i = eqline( A, B, m1, d1 );

    m2i = eqline( C, D, m2, d2 );

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

	/*
		空间直线的两点式： 
		（类似于平面坐标系中的两点式） 
		设两点为A(x1,y1,z1),B（x2,y2,z2)
		则直线AB方程为(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
	 */
	if (A[2] - B[2] == 0)
	{
		I[2] = A[2];
	}
	else if (A[0] - B[0] != 0)
	{
		//(x-x1)/(x2-x1)=(z-z1)/(z2-z1)
		I[2] = (A[2] - B[2])*(I[0]-B[0])/(A[0] - B[0]) + B[2];
	}
	else
	{
		// (y-y1)/(y2-y1)=(z-z1)/(z2-z1)
		I[2] = (A[2] - B[2])*(I[1]-B[1])/(A[1] - B[1]) + B[2];
	}
#endif

}


#define INSERT_SEGMENTS 10.0
/* 求取空间两个点A和B之间在直线AB上的的若干(10)个间隔相等的点，保存在vec_list里面  */
void interpolation(osg::Vec3 A, osg::Vec3 B, std::vector<osg::Vec3> &vec_list, float segs)
{
	double delx, dely, delz, step, delta;
	int count;

	delx = B[0]-A[0];
	dely = B[1]-A[1];
	delz = B[2]-A[2];
	step = segs;

	if (fabs(delx)>fabs(dely))
	{
		count = fabs(delx) / step + 1;
		delta = delx / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3 iVec;
			iVec[0] = A[0] + delta*(i+1);
			iVec[1] = (delta*(i+1))*dely/delx + A[1];
			iVec[2] = (delta*(i+1))*delz/delx + A[2];
			vec_list.push_back(iVec);
		}
	}
	else if (fabs(dely)>fabs(delx))
	{
		count = fabs(dely) / step + 1;
		delta = dely / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3 iVec;
			iVec[1] = A[1] + delta*(i+1);
			iVec[0] = (delta*(i+1))*delx/dely + A[0];
			iVec[2] = (delta*(i+1))*delz/dely + A[2];
			vec_list.push_back(iVec);
		}
	}
	else if((fabs(delx) == 0.0) && (fabs(dely) == 0.0))
	{
		osg::Vec3 iVec;
		iVec[0] = A[0];
		iVec[1] = A[1];
		iVec[2] = A[2];
		vec_list.push_back(iVec);
	}
	else
	{
		/* fabs(delx)和fabs(dely)相等，但都不等于0.0 */
		count = fabs(delx) / step + 1;
		delta = delx / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3 iVec;
			iVec[0] = A[0] + delta*(i+1);
			iVec[1] = (delta*(i+1))*dely/delx + A[1];
			iVec[2] = (delta*(i+1))*delz/delx + A[2];
			vec_list.push_back(iVec);
		}
	}
}


double LenOf2V(const osg::Vec3 v0, const osg::Vec3 v1)
{ return sqrt( (v0.x() - v1.x()) * (v0.x() - v1.x()) + (v0.y() - v1.y()) * (v0.y() - v1.y()) + (v0.z() - v1.z()) * (v0.z() - v1.z())); }

double LenOf2VF(const osg::Vec3 v0, const osg::Vec3 v1)
{ return sqrt( (v0.x() - v1.x()) * (v0.x() - v1.x()) + (v0.y() - v1.y()) * (v0.y() - v1.y())); }

/* 找出所有的奇点，奇点的定义是：在一个线环中出现的次数不等于2的点，一个纯线环上的每个点应该都只出现2次。 */
void FindOutOddPoints(std::vector<SimEdge> *edges)
{
	/* 这一段程序为了删除线环里面多出来的套圈，使得线环变成一个纯的线环。*/
	std::vector<VerIdx> AllPoints;
	std::vector<VerIdx> OddPoints;
	std::vector<VerPair> OddPointPairs;

	/* 首先，如果有重复的边，先去掉。*/
	std::vector<SimEdge> remain;
	remain.push_back(edges->front());
	std::vector<SimEdge>::iterator p = edges->begin() + 1;
	for( ; p != edges->end(); p++)
	{
		int iFlag = 0;
		for(std::vector<SimEdge>::iterator q = remain.begin(); q != remain.end(); q++)
		{
			if(*p == *q)	{	iFlag = 1;	break;	}
		}
		if(!iFlag)	remain.push_back(*p);
	}
	edges->clear();
	for(std::vector<SimEdge>::iterator q = remain.begin(); q != remain.end(); q++)	edges->push_back(*q);
		
	/* 第一步，找出边集合中所有的定点 */
	while(1)
	{
		int index = 0;
		AllPoints.clear();
		for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
		{
			int iFlag = 0;
			osg::Vec3 Pt = p->A;
			for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
			{
				if(r->Ver == Pt)	{	r->times++; iFlag = 1; break; }
			}
			if(!iFlag)
			{
				VerIdx vi; vi.Ver = Pt; vi.idx = index++; vi.times = 1;
				AllPoints.push_back(vi);
			}
			
			iFlag = 0; Pt = p->B;
			for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
			{
				if(r->Ver == Pt)	{	r->times++; iFlag = 1; break; }
			}
			if(!iFlag)
			{
				VerIdx vi; vi.Ver = Pt; vi.idx = index++; vi.times = 1;
				AllPoints.push_back(vi);
			}
		}
		
		/* 第二步要找出所有的奇点，奇点的定义是：在一个线环中出现的次数不等于2的点，一个纯线环上的每个点应该都只出现2次。 */
		OddPoints.clear();
		for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
		{
			if(r->times == 1)
			{
				VerIdx vi; vi.Ver = r->Ver; vi.idx = r->idx; vi.times = r->times;
				OddPoints.push_back(vi);
			}
		}
		
		/* 增加一例特殊处理状况：如果有且只有两个OddPoints，说明这个线环是非封闭的线环。*/	
		/* 把这两个奇点连接起来，使得线环变成一个封闭线环。2012-10-27 */
		if(OddPoints.size() == 2)
		{
			/* 把这两个奇点连接成一个段，*/
			SimEdge se;	 
			se.A = OddPoints[0].Ver; se.B = OddPoints[1].Ver; se.x = 0; se.y = 0; se.z = 0;
			edges->push_back(se);
			OddPoints.clear();
		}

		/* 看看奇点之间有没有直接连接的线段，有就把这个线段去掉。 */
		if(OddPoints.size() > 0)
		{
			std::vector<SimEdge> errorEdges;
			for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
			{
				for(unsigned int j=0; j<OddPoints.size(); j++)
				{
					if((OddPoints[j].times == 1) && (( p->A == OddPoints[j].Ver) || ( p->B == OddPoints[j].Ver)))
					{
						SimEdge se;	 
						se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
						errorEdges.push_back(se);
					}
				}
			}
			
			/* 重新整理totalEdges */
			std::vector<SimEdge> temp;
			for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
			{
				int iFlag = 0;
				for( std::vector<SimEdge>::iterator q = errorEdges.begin(); q != errorEdges.end(); q++ )
				{
					if((*q == *p))	{	iFlag = 1; break; }	
				}
				if(!iFlag) 
				{
					SimEdge se;	 
					se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
					temp.push_back(se);
				}
			}
			
			edges->clear();
			for( std::vector<SimEdge>::iterator p = temp.begin(); p != temp.end(); p++ )
			{
				SimEdge se;	 
				se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
				edges->push_back(se);
			}
		}else break;
	}	//while
}


/* 找出断头点并逐步消除断头。 */
void RemoveOddPoints(std::vector<SimEdge> *edges)
{
	/* 这一段程序为了删除线环里面多出来的套圈，使得线环变成一个纯的线环。*/
	std::vector<VerIdx> AllPoints;
	std::vector<VerIdx> OddPoints;
	std::vector<VerPair> OddPointPairs;
		
	/* 第一步，找出边集合中所有的定点 */
	while(1)
	{
		int index = 0;
		AllPoints.clear();
		for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
		{
			int iFlag = 0;
			osg::Vec3 Pt = p->A;
			for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
			{
				if(r->Ver == Pt)	{	r->times++; iFlag = 1; break; }
			}
			if(!iFlag)
			{
				VerIdx vi; vi.Ver = Pt; vi.idx = index++; vi.times = 1;
				AllPoints.push_back(vi);
			}
			
			iFlag = 0; Pt = p->B;
			for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
			{
				if(r->Ver == Pt)	{	r->times++; iFlag = 1; break; }
			}
			if(!iFlag)
			{
				VerIdx vi; vi.Ver = Pt; vi.idx = index++; vi.times = 1;
				AllPoints.push_back(vi);
			}
		}
		
		/* 第二步要找出所有的奇点，奇点的定义是：在一个线环中出现的次数不等于2的点，一个纯线环上的每个点应该都只出现2次。 */
		OddPoints.clear();
		for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
		{
			if(r->times == 1)
			{
				VerIdx vi; vi.Ver = r->Ver; vi.idx = r->idx; vi.times = r->times;
				OddPoints.push_back(vi);
			}
		}
		
		/* 看看奇点之间有没有直接连接的线段，有就把这个线段去掉。 */
		if(OddPoints.size() > 0)
		{
			std::vector<SimEdge> errorEdges;
			for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
			{
				for(unsigned int j=0; j<OddPoints.size(); j++)
				{
					if((OddPoints[j].times == 1) && (( p->A == OddPoints[j].Ver) || ( p->B == OddPoints[j].Ver)))
					{
						SimEdge se;	 
						se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
						errorEdges.push_back(se);
					}
				}
			}
			
			/* 重新整理totalEdges */
			std::vector<SimEdge> temp;
			for( std::vector<SimEdge>::iterator p = edges->begin(); p != edges->end(); p++ )
			{
				int iFlag = 0;
				for( std::vector<SimEdge>::iterator q = errorEdges.begin(); q != errorEdges.end(); q++ )
				{
					if((*q == *p))	{	iFlag = 1; break; }	
				}
				if(!iFlag) 
				{
					SimEdge se;	 
					se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
					temp.push_back(se);
				}
			}
			
			edges->clear();
			for( std::vector<SimEdge>::iterator p = temp.begin(); p != temp.end(); p++ )
			{
				SimEdge se;	 
				se.A = p->A; se.B = p->B; se.x = p->x; se.y = p->y; se.z = p->z;
				edges->push_back(se);
			}
		}else break;
	}	//while
}



/* 把 validBorderSegment 组成一个线环 */
void CreateLineLoopFromValidBorderSegment(osg::ref_ptr<osg::Vec3Array> coords)
{
	std::vector<SimEdge>  totalEdges;	totalEdges.clear();
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		SimEdge se;	 
		se.A = q->pS; se.B = q->pE; se.x = 0; se.y = 0; se.z = 0; se.keep = true;
		if(se.A != se.B)	totalEdges.push_back(se);
	}

	/* 找到并消除奇点对*/
	FindOutOddPoints(&totalEdges);
	
	/* 有可能 validBorderSegment 集合不是一个封闭的段，那么这时候 totalEdges 的数目就会变成0 */
	if(totalEdges.size() == 0)	{		return;	}

	while(totalEdges.size() > 2)
	{
		/* 把所有点连成一条线LINE_LOOP*/
		for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ ) p->keep = true;
		SimEdge start = totalEdges.front();

CPA_LOOP_0:
		SimEdge e = start;
		osg::Vec3 A = e.A;
		osg::Vec3 B = e.B;
		osg::Vec3 firstA = A;
		bool bEnd = true;	//如果当前的所有边凑不成一个线环，bEnd就会等于true，导致最终退出do..while循环。如果能凑成一个线环，bEnd始终为false.
		do {
			int i = 0;
			bEnd = true;
			for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
			{
				if( p->has(B) && !(e == *p) && p->keep )
				{
					coords->push_back( A );
					A = B;
					if( B == p->A )
						B = p->B;
					else
						B = p->A;
	
					p->keep = false;
					e = *p;
					bEnd = false;
					break;
				}
			}
			if (bEnd)
			{
				if((start.has(e.A)) || (start.has(e.B)))
				{	
					if(start.has(e.B))	coords->push_back( A );
					e = start;	bEnd = false;
				}
			}
		}while( !(e == start) && !bEnd );
	
		/* 如果上面循环是以bEnd为条件退出循环体的话，说明得到的这个线环并不是合格的线环。2012-10-11 */
		if(e == start)	coords->push_back(firstA);
		else
		{
			/* 把 e 这个边，从原有totalEdges 里面拿掉，然后重新计算。 */
			std::vector<SimEdge> remain;
			for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
			{
				if (!(*p == e))	remain.push_back(*p);
			}
			totalEdges.clear();
			for( std::vector<SimEdge>::iterator p = remain.begin(); p != remain.end(); p++ )
			{
				p->keep = true;
				totalEdges.push_back( *p );
			}


			/* 重新找LINELOOP */
			coords->clear();
			goto CPA_LOOP_0;
		}
		
		/* 找到一个线环之后，先把这个线环存起来，存到 multiDc*/
		if(coords->size() > 3)	
		{
			osg::ref_ptr<osg::Vec3Array> nCords = dynamic_cast<osg::Vec3Array *> (coords->clone(osg::CopyOp::DEEP_COPY_ALL));
			multiDc.push_back(nCords);
		}
		coords->clear();

		/* 把 totalEdges 里面 keep == false的边移走一个，然后重新计算。 */
		{
			/* 先移走一个false边 */
			std::vector<SimEdge>::iterator remove;
			int iFalse = 0;
			for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
			{
				if (p->keep == false)	{	remove = p;	iFalse++;}
			}
			if(iFalse > ( totalEdges.size() - 3) ) break;

			std::vector<SimEdge> remain;
			for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
			{
				if (p != remove) remain.push_back(*p);
			}
			/* 剩下的边先去掉 奇点 */
			RemoveOddPoints(&remain);
			/* 剩下的边再交还给totalEdges 重新计算下一个线环*/
			totalEdges.clear();
			if(remain.size() > 3)
			{
				for( std::vector<SimEdge>::iterator p = remain.begin(); p != remain.end(); p++ )
				{
					totalEdges.push_back( *p );
				}
			}else break;
		}
	}	// while(totalEdges.size() > 0)
}

typedef struct _SEGMENTINTERSECTION_
{
	int index;		//与文化信息相交的边界段的序号；
	int num;		//文化信息与这段边界相交了几次；
}SegmentIntersection;

typedef struct _DISTANCEOFSEGMENT_
{
	osg::Vec3	vertex;		//交点
	double		distance;	//交点与段端点的距离。
}DistanceOfSegment;

/* 这里要处理一种情况，如果文化信息曲线与一段边界有多次的相交，该怎么办？*/
std::vector<int> removeValidSegmentWhenWithMultiIntersections;
std::vector<FindValidVertex> intsResult;  		//记录了文化信息外围环与地形边界相交的点集合。
void CheckSegmentsIntersection()
{
	unsigned int i,j,k,m;
	std::vector<SegmentIntersection> allSis;	allSis.clear();
	
	/* 首要检查一下 removeValidSegmentWhenWithMultiIntersections 里面是否真的有多条线段相交。*/
	/* 因为判断removeValidSegmentWhenWithMultiIntersections的条件是：*/
	/* 边界：n-1点相交，n点不相交，说明 n-1 和 n 两点之间可能有线段通过，这仅仅是可能。*/
	/* 此时如果 n 点与文化信息交点重合，那么就会多增加一次n-1到n的线段，这样会与可能的原有相交边界线段部分重合 */
	/* 但是不能保证 n-1 和 n 两点之间确实有线段通过，所以还要检查一下。*/
	if(removeValidSegmentWhenWithMultiIntersections.size() > 0)
	{
		std::vector<int> tmpMI;
		for(unsigned int i=0; i<removeValidSegmentWhenWithMultiIntersections.size(); i++)
		{
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
			{
				if(q->idxE == removeValidSegmentWhenWithMultiIntersections[i])
				{
					tmpMI.push_back(removeValidSegmentWhenWithMultiIntersections[i]);
					break;
				}
			}
		}
		removeValidSegmentWhenWithMultiIntersections.clear();
		if(tmpMI.size() > 0)
		{
			for(unsigned int i=0; i<tmpMI.size(); i++)	removeValidSegmentWhenWithMultiIntersections.push_back(tmpMI[i]);
		}
	}	

	/* 首先在相交集合 validBorderSegment 里面统计区段 */
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		if(allSis.size() > 0)
		{
			int iFlag = 0;
			for( std::vector<SegmentIntersection>::iterator r= allSis.begin(); r != allSis.end(); r++ )
			{
				if(q->border == r->index)
				{	iFlag = 1; r->num++;	break;}
			}
			if(!iFlag)
			{	SegmentIntersection si; si.index = q->border; si.num = 1; allSis.push_back(si);	}
		}else{
			SegmentIntersection si; si.index = q->border; si.num = 1; allSis.push_back(si);
		}
	}
	
	/* 判断统计结果allSis的size()与validBorderSegment的size()是否相等，
	   如果相等，说明每一个相交段只有一个交点，是正常情况不作处理。
	   如果不相等，说明某些相交段内有多个交点。
	 */
	if( allSis.size() < validBorderSegment.size() )
	{
		for(i=0; i<allSis.size(); i++)
		{
			if(allSis[i].num > 1)
			{
				/* 找出多点相交段的所有相交点，并计算出该交点和段端点的距离，收集到集合allDOS里面 */
				std::vector<LineSeg> currBorderSegment; currBorderSegment.clear();
				std::vector<DistanceOfSegment> allDOS;	allDOS.clear();
				osg::Vec3 S_fst, S_next, S_end;	S_fst.set(0.0, 0.0, 0.0);
				int IntFlag;
				for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				{
					if(q->border == allSis[i].index)
					{
						if((S_fst.x() == 0.0) && (S_fst.y() == 0.0))
						{
							/* 找出这一边界段的原始的两个端点，保存在S_fst和S_end里面 */
							S_fst = q->pE;
							if(borderIntersectLS[q->border].pS == q->pE) S_end = borderIntersectLS[q->border].pE; else S_end = borderIntersectLS[q->border].pS;
							IntFlag = q->validE;
						}
						
						DistanceOfSegment dos;	dos.vertex = q->pS; dos.distance = LenOf2V(q->pS, S_fst);	allDOS.push_back(dos);
						
					}else{
						LineSeg ls;
						ls.pS = q->pS;   ls.idxS = q->idxS;	ls.validS = q->validS;
						ls.pE = q->pE;   ls.idxE = q->idxE;	ls.validE = q->validE;
						ls.border = q->border;
						currBorderSegment.push_back(ls);						
					}
				}
				
				/* IntFlag == 2 表示原始边界段两个端点都与文化信息相交 */
				/* 还要判断这个边界段两个端点是否相交状况都为1，如果都为1，记录下这个段，以备后边把它去掉。*/		
				if(IntFlag == 2)
				{
					DistanceOfSegment dos;	dos.vertex = S_end; dos.distance = LenOf2V(S_fst, S_end);	allDOS.push_back(dos);
					removeValidSegmentWhenWithMultiIntersections.push_back(allSis[i].index);
				}

				/* 对这个集合里面的所有元素，按distance从小到大进行排序 */
				for(k=0; k<(int)(allDOS.size()-1); k++)
				{
					for(j=k+1; j<(int)(allDOS.size()); j++)
					{
						if(allDOS[k].distance > allDOS[j].distance)
						{
							DistanceOfSegment tmp;
							tmp.vertex = allDOS[k].vertex; tmp.distance = allDOS[k].distance;
							allDOS[k].vertex = allDOS[j].vertex; allDOS[k].distance = allDOS[j].distance;
							allDOS[j].vertex = tmp.vertex; allDOS[j].distance = tmp.distance;
						}
					}
				}
				
				/* 排完序之后，将有效段加入到 validBorderSegment 里面 
				   分三步，第一步，先把计算出来的有效段加入到临时集合currBorderSegment里面
				   第二步，把原来 validBorderSegment 里面的无效段全部删除；
				   第三步，把 validBorderSegment 和 currBorderSegment 有效的组合在一起 */

				m = 0;
				unsigned int mod;
				if(IntFlag > 0) mod = 1; else mod = 0;
				while(m < allDOS.size())
				{
					S_next = allDOS[m].vertex;
					if(mod)
					{
						LineSeg ls;
						ls.pS = S_fst;   ls.idxS = -1;	ls.validS = 1;
						ls.pE = S_next;  ls.idxE = allSis[i].index;	ls.validE = 1;
						ls.border = allSis[i].index;
						currBorderSegment.push_back(ls);
					}
					S_fst = S_next; mod = (mod+1)%2; m++;
				}
				
				/* 整合 */
				validBorderSegment.clear();
				for( std::vector<LineSeg>::iterator q= currBorderSegment.begin(); q != currBorderSegment.end(); q++ )
				{
					LineSeg ls;
					ls.pS = q->pS;   ls.idxS = q->idxS;	ls.validS = q->validS;
					ls.pE = q->pE;   ls.idxE = q->idxE;	ls.validE = q->validE;
					ls.border = q->border;
					validBorderSegment.push_back(ls);						
				}		

			}
		}
	}

	/* 如果一个有效线段的两个端点相同，去掉这个线段。*/
	std::vector<LineSeg> tmpSegment;
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		double dist = LenOf2V(q->pE, q->pS);
		if(dist > MINDISTTOMERGEVERTEX) tmpSegment.push_back(*q);
	}
	validBorderSegment.clear();
	for( std::vector<LineSeg>::iterator q= tmpSegment.begin(); q != tmpSegment.end(); q++ )
		validBorderSegment.push_back(*q);

}

/* 检查特殊情况。 */
typedef struct _BorderVertexStatus_
{
	int idx;
	int isValid;	
}BorderVertexStatus;
std::vector<int> mergeVertex;			//保存边界段端点和交点的合并状况

void CheckSpecialStatus()
{
	/* Add by Galen 2013-02-05 */
	/* 如果 validBorderSegment 内有两段相同，说明有两个文化信息段与地形边的交点重合。这个交点虽然落在了地形边上，但是并未对这个地形边进行实质性的切割。*/
	/* 应当把这两个相同段去掉。并且把那个交点也从交点集合 allBorderVertices 里面去掉，*/
	int hasSameSegment = 0;
	std::vector<osg::Vec3> interpVFs;
	for( std::vector<LineSeg>::iterator q = validBorderSegment.begin(); q != validBorderSegment.end(); ++q )
	{
		/* 先复制一份 */
		std::vector<LineSeg> tmpSegment;
		for( std::vector<LineSeg>::iterator p = validBorderSegment.begin(); p != validBorderSegment.end(); ++p )
			if(p != q)	tmpSegment.push_back(*p);
		
		/* 再比较, 注意：！！！这里应当再多加一个判断条件，就是：交点所在边段两个端点与文化信息面的相交状态应该是相同的。稍后再加 2013-02-06 */
		for( std::vector<LineSeg>::iterator p = tmpSegment.begin(); p != tmpSegment.end(); ++p )
		{
			if((p->pS == q->pS) && (p->pE == q->pE))
			{
				q->border = -2;	hasSameSegment = 1;	interpVFs.push_back(p->pS);
			}
		}
	}
	if(hasSameSegment)
	{
		/* 先把相同边从validBorderSegment里面去掉 */
		std::vector<LineSeg> tmpSegment;
		for( std::vector<LineSeg>::iterator p = validBorderSegment.begin(); p != validBorderSegment.end(); ++p )
			tmpSegment.push_back(*p);
		validBorderSegment.clear();
		for( std::vector<LineSeg>::iterator p = tmpSegment.begin(); p != tmpSegment.end(); ++p )
		{
			if(p->border != -2)
			{
				validBorderSegment.push_back(*p);
			}
		}

		/* 再从交点集合里面去掉相关的交点。*/
		std::vector<FindValidVertex> tmpBorderVertices;	
		std::copy(allBorderVertices.begin(), allBorderVertices.end(), back_inserter(tmpBorderVertices));
		allBorderVertices.clear();
		for( std::vector<FindValidVertex>::iterator p = tmpBorderVertices.begin(); p != tmpBorderVertices.end(); ++p )
		{
			int iFlag = 0;
			for(std::vector<osg::Vec3>::iterator v = interpVFs.begin(); v != interpVFs.end(); ++v)	
			{
				if(p->pVF == *v) {	iFlag = 1;	break;	}
			}
			if(!iFlag)	allBorderVertices.push_back(*p);
		}
	}
	/* Add End */

	/* 修正一下validBorderSegment里面的 validE 值。仅在下列状况下：*/
	/* validBorderSegment[].validE = 0(先计算出来的)，但这个 pE点也是交点(后计算出来的)，那么要把这个 validBorderSegment[].validE 改为 1 */
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		if(q->validE == 0)
		{
			int curIdx = q->idxE;
			if(terrainBorderVertices[curIdx].newIdx == 1)
				q->validE = 1;
		}
	}

	/* 特殊情况：如果只有两个交点，但在不同的段内，而且这两个不同的段分属两个不同的边，那么这两个交点一定夹着一个角。
	   改正办法：把原本不相连接的两个段连起来 */
	if(	(validBorderSegment.size() == 2) && 
		(validBorderSegment[0].border != validBorderSegment[1].border) &&
		(validBorderSegment[0].validE == 0) &&
		(validBorderSegment[1].validE == 0) )
	{
		osg::Vec3 connVer;
		for(int ic=0;ic<4;ic++)
		{
			if(conner[ic].pVF == validBorderSegment[0].pE) {	validBorderSegment[1].pE = validBorderSegment[0].pE;	break;	}
			if(conner[ic].pVF == validBorderSegment[1].pE) {	validBorderSegment[0].pE = validBorderSegment[1].pE;	break;	}
		}
		return;
	}


	/* 特殊情况：如果只有两个交点，还都在一个段内，那么，这个段的首尾两个端点，相交状态应当是一致的，如果不是一致的就错了。
	   改正办法：如果validE状态是1，那么把validE状态改成0. */
	if(	(validBorderSegment.size() == 2) && 
		(validBorderSegment[0].border == validBorderSegment[1].border) &&
		(validBorderSegment[0].validE == 1) &&
		(validBorderSegment[1].validE == 1) )
	{
		validBorderSegment[0].validE = validBorderSegment[1].validE = 0;
	}

	/* 特殊情况：求出来的有效相交段个数为0，但是半相交段的个数却不为0，说明相交面积太小了，可忽略不计。
	   改正办法：调整边界的各个端点的相交状态。
	   修正：如果相交点在不同的边，上述的特殊情况就不存在。另外有效线段为0不等于相交点个数为0，这里有都和边界点重合的情况。
	   只有两个相交点在同一边上，上述特殊情况才成立。*/
	if((borderHalfIntersectLS.size() > 0) && ( validBorderSegment.size() == 0))
	{
		/* 先判断交点情况，intsResult[]，如果至少有两个交点而且这两个交点同不同的边相交，那么上述特殊情况就不成立。*/
		if(!((intsResult.size() >= 2) && (intsResult[0].isValid != intsResult[1].isValid)))
		{
			/* 当前边界点中有多少个相交点*/
			int intBordVers = 0;
			for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
			{
				if(q->newIdx != 0) intBordVers++;
			}
	
			/* 如果当前边界点中相交点个数少于 10 个。*/
			if(intBordVers < 10)
			{
				/* 把边界点全设成相交点*/
				for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
				{
					q->newIdx = 0;
				}
			}
	
			/* 清除掉所有的有效相交线段 */	
			validBorderSegment.clear();
			borderHalfIntersectLS.clear();
			return;
		}
	}

	/* 特殊情况：如果相交段的个数大于交点的个数，说明有交点计算不出来，文化信息一条边和边界近似重合。
	   改正办法：调整边界的各个端点的相交状态。 */
	if(borderHalfIntersectLS.size() >  (validBorderSegment.size() + mergeVertex.size()))
	{
		/* 如果剩余的两个半相交段相邻，且分属不同的边，说明这两个半相交段所夹的是一个角点。*/
		/* 这个角点的两边相交段没有计算出交点，说明这个角点正巧落在了文化信息线上。那么把这个角点还原即可。*/
		std::vector<LineSeg> validHalfIntersectLS;
		for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
		{
			if(r->border != -2)
			{
				validHalfIntersectLS.push_back(*r);
			}
		}
		/* 只剩下两个有效的半相交段 */
		if(validHalfIntersectLS.size() == 2)
		{
			/* 这两个有效半相交段要相邻 */
			if ((validHalfIntersectLS[0].idxS == validHalfIntersectLS[1].idxE) ||
					(validHalfIntersectLS[0].idxE == validHalfIntersectLS[1].idxS) )
			{
				/* 这两个有效半相交段分属不同的边 */
				if(validHalfIntersectLS[0].border != validHalfIntersectLS[1].border)
				{
					if(validHalfIntersectLS[0].idxS == validHalfIntersectLS[1].idxE)
					{
						int idx = validHalfIntersectLS[0].idxS;
						terrainBorderVertices[idx].newIdx = 1 - terrainBorderVertices[idx].newIdx;
					}
					if(validHalfIntersectLS[0].idxE == validHalfIntersectLS[1].idxS)
					{
						int idx = validHalfIntersectLS[1].idxS;
						terrainBorderVertices[idx].newIdx = 1 - terrainBorderVertices[idx].newIdx;
					}
					return;
				}
			}
		}
		
		/* 先看看是否所有内部点都和文化信息面相交，如果相交，就说明有重合现象。*/
		if(AllInnerVertexIntersection)
		{
			/* 当前边界点中有多少个相交点*/
			int intBordVers = 0;
			for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
			{
				if(q->newIdx != 1) intBordVers++;
			}

			/* 如果当前边界点中不相交点个数少于 3 个。*/
			if(intBordVers < 3)
			{
				/* 把边界点全设成相交点*/
				for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
					q->newIdx = 1;
				/* 清除掉所有的有效相交线段 */	
				validBorderSegment.clear();
			}
			return;
		}
		
		/* 对于有一部分内部点与文化信息不相交的情况 */
		std::vector<int> toModifySegment;
		for(unsigned int i=0; i<borderHalfIntersectLS.size(); i++)
		{
			int iFlag = 0;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
			{
				if((q->idxE == borderHalfIntersectLS[i].idxS) || (q->idxE == borderHalfIntersectLS[i].idxE))
				{ iFlag = 1; break; }
			}
			if(!iFlag)	toModifySegment.push_back(i);
		}
		
		if(toModifySegment.size())
		{
			for(unsigned int i=0; i<toModifySegment.size(); i++)
			{
				int iFlag = 0;
				for(unsigned int j=0; j<mergeVertex.size(); j++)
				{
					int idx = borderHalfIntersectLS[toModifySegment[i]].idxS;
					if (( abs(idx - mergeVertex[j]) <= 1) && (borderHalfIntersectLS[toModifySegment[i]].border != -2))
					{
						iFlag = 1; break;
					}
				}
				
				if(!iFlag)
				{
					int idx = borderHalfIntersectLS[toModifySegment[i]].idxS;
					int num = terrainBorderVertices.size();
					if(terrainBorderVertices[idx].newIdx == 0)
					{
						terrainBorderVertices[idx].newIdx = 1;
					}

					idx = borderHalfIntersectLS[toModifySegment[i]].idxE;
					if(terrainBorderVertices[idx].newIdx == 0)
					{
						terrainBorderVertices[idx].newIdx = 1;
					}
				}
			}
			
			/* 去掉无效的validBorderSegment边 */
			std::vector<LineSeg> temp;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
			{
				if(q->validE != 2)
				{ 
					LineSeg ls;
					ls.pS = q->pS;   ls.idxS = q->idxS;	ls.validS = q->validS;
					ls.pE = q->pE;   ls.idxE = q->idxE;	ls.validE = q->validE;
					ls.border = q->border;
					temp.push_back(ls);
				}
			}
			
			validBorderSegment.clear();
			for( std::vector<LineSeg>::iterator q= temp.begin(); q != temp.end(); q++ )
			{
				LineSeg ls;
				ls.pS = q->pS;   ls.idxS = q->idxS;	ls.validS = q->validS;
				ls.pE = q->pE;   ls.idxE = q->idxE;	ls.validE = q->validE;
				ls.border = q->border;
				validBorderSegment.push_back(ls);
			}
		}
		return;
	}

	/* 特殊情况：如果一个边界段的两个端点都与文化信息相交或者不相交，但是这个边界段内却有而只有一个交点，这种情况显然是错误的。但是却出现了!!!.
	   改正办法：找到最近的相交段，修改本段和最近相交段之间的相交信息。 */

	/* 首先要找到两个端点相交状态相同，但是段内却有交点的所有段 */
	std::vector<int> bothSegment;
	std::vector<int> diffSegment;
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		if(q->validE != 1) bothSegment.push_back(q->border);
		if(q->validE == 1) diffSegment.push_back(q->border);
	}
	
	/* 处理边界段两个端点是不同相交状态的特殊情况 */
	if(diffSegment.size())
	{
		/* 其次要找到两个端点相交状态不同，但是段内却有偶数个交点的所有错误段 */
		std::vector<int> errorSegment;
		for(unsigned int i=0; i<diffSegment.size(); i++)
		{
			int iNum = 0;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				if(q->border == diffSegment[i]) iNum++;
			if(iNum%2 == 0) errorSegment.push_back(diffSegment[i]);
		}
		
		/* 然后对这些错误段进行处理 */
		if(errorSegment.size())
		{
			/* 根据errorSegment判断相交状态 */
			for(unsigned int i=0; i<errorSegment.size(); i++)
			{
				for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				{
					if(q->border == errorSegment[i])
					{
						/* 调整validBorderSegment 相交状态 */
						/* 调整terrainBorderVertices对应值的 相交状态 */
						if(q->idxE != q->border)
						{
							if(terrainBorderVertices[q->border].newIdx == 0)	
							{	
								q->validE = 0;
								terrainBorderVertices[q->border + 1].newIdx = 0;
							} else 	{
								terrainBorderVertices[q->border + 1].newIdx = 1;
								q->validE = 2;
							}
						}else{
							if(terrainBorderVertices[q->border + 1].newIdx == 0)
							{
								q->validE = 0;
								terrainBorderVertices[q->border].newIdx = 0;
							} else 	{
								q->validE = 2;
								terrainBorderVertices[q->border].newIdx = 1;
							}
						}
					}
				}
			}
		}		
	}
	
	/* 处理边界段两个端点是相同相交状态的特殊情况 */
	if(bothSegment.size())
	{
		/* 其次要找到两个端点相交状态相同，但是段内却有奇数个交点的所有错误段, 增加一个限定条件：这个奇数值大于1。  */
		std::vector<int> errorSegment;
		for(unsigned int i=0; i<bothSegment.size(); i++)
		{
			int iNum = 0;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				if(q->border == bothSegment[i]) iNum++;
			if((iNum%2 == 1)&&(iNum > 1)) errorSegment.push_back(bothSegment[i]);
		}

		/* 然后对这些错误段进行处理 */
		if(errorSegment.size())
		{
			/* 判断所有半相交边界段，是否有没求出交点的半相交边界段。*/

			/* 找出所有半相交段里面没有交点的段的两个端点及状态值。*/
			std::vector<BorderVertexStatus> allBVS; allBVS.clear();
			for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
			{
				if(r->border != -2)
				{
					BorderVertexStatus bvs;
					bvs.idx = r->idxS ; bvs.isValid = r->validS; allBVS.push_back(bvs);
					bvs.idx = r->idxE ; bvs.isValid = r->validE; allBVS.push_back(bvs);
				}				
			}

			/* 根据errorSegment判断相交状态 */
			for(unsigned int i=0; i<errorSegment.size(); i++)
			{
				int currFlag = terrainBorderVertices[errorSegment[i]].newIdx;
				int currIdx  = terrainBorderVertices[errorSegment[i]].idx;
				int minDist  = 999, minIdx = -1;
				
				/* 找到与当前错误相交段最短距离的半相交边界段的端点 */
				for(unsigned int j=0; j<allBVS.size(); j++)
				{
					if(currFlag == allBVS[j].isValid)
					{
						int dist = abs(currIdx - allBVS[j].idx);
						if(dist < minDist) {	minDist = dist; minIdx = allBVS[j].idx ;	}
					}
				}
				
				/* 转换terrainBorderVertices里面部分点的相交状态。 */
				unsigned int m, n;
				if(minIdx <= currIdx) {	m = minIdx;	n = currIdx;	}
				if(minIdx > currIdx) {	m = currIdx + 1;	n = minIdx;	}
				for(unsigned int j=m; j<=n; j++)
				{
					if(currFlag==0)	terrainBorderVertices[j].newIdx = 1;
					if(currFlag==1)	terrainBorderVertices[j].newIdx = 0;
				}
				
				/* 最后改变validBorderSegment里面相应段的状态 */
				for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				{
					if((q->validE != 1) && (q->border == errorSegment[i]))
					{
						q->validE = 1;	
						if(terrainBorderVertices[currIdx].newIdx == 1) { q->pE = terrainBorderVertices[currIdx].pVF; q->idxE = currIdx;	}	else  { q->pE = terrainBorderVertices[currIdx+1].pVF; q->idxE = currIdx + 1;	}	
						break;
					}
				}
			}
		}
	}
}

/* 如果有一个交点落在了角点上 */
void CheckSpecialIntVertex(osg::Vec3 I)
{
	int iFlag = 0;
	for(unsigned int i=0; i<4; i++)
	{
		double dist = LenOf2V(I, conner[i].pVF);
		if(dist < MINDISTTOMERGEVERTEX)	{ iFlag = 1; break;}
	}
	if(!iFlag) return;
	
	/* 找出这个角点所在的半相交段。*/
	std::vector<LineSeg> validHalfIntersectLS;
	for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
	{
		double distS = LenOf2V(I, r->pS);
		double distE = LenOf2V(I, r->pE);
		if((distS < MINDISTTOMERGEVERTEX)||(distE < MINDISTTOMERGEVERTEX))
		{
			validHalfIntersectLS.push_back(*r);
		}
	}

	/* 如果不是相邻两个半相交段都包含这个角点，则退出。 */
	if(validHalfIntersectLS.size() != 2) return;
#if 0
	/* 这段代码的意图不知何故。忘记了，当初为什么要把“经过角点的相交段”给删掉。先关掉，2012-10-19*/
	std::vector<LineSeg> tmpSegment;
	for(unsigned int m = 0; m<validBorderSegment.size(); m++)
	{
		int jFlag = 0;
		for(std::vector<LineSeg>::iterator r = validHalfIntersectLS.begin(); r != validHalfIntersectLS.end(); r++)
		{
			if((validBorderSegment[m].idxE == r->idxS) || (validBorderSegment[m].idxE == r->idxE))
			{
				jFlag = 1; break;
			}
		}
		if(!jFlag) tmpSegment.push_back(validBorderSegment[m]);
	}
	validBorderSegment.clear();
	for(unsigned int m = 0; m<tmpSegment.size(); m++)
		validBorderSegment.push_back(tmpSegment[m]);
#endif
}


/* 在 checkBorderVertices 里面查找一个顶点，满足：A:onCulture = 1; B:这个点和输入参数点距离小于1.0 .*/
/* 如果找到这样的点，把它的onCulture = 0 */
void FindVertexInCBVsets(osg::Vec3 vx)
{
	for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
	{
		if((p->onDetailed == 0) && (p->onCulture == 1)) 
		{
			double dist = LenOf2V(vx, p->fvv.pVF);
			if(dist < MINDISTTOMERGEVERTEX)
			{
				p->onCulture = 0; break;
			}
		}
	}
}

/* 当有合并点的情况出现后，要把 faceVertice 集合里面的相应的点要修改一下。否则，后面做德罗尼三角化的时候就会得不到正确结果。*/
void ModifyFaceVertexWhenMerge(FindValidVertex fvv)
{
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
	{
		double dist = LenOf2VF(fvv.pVF, f->pVF);
		if(dist < MINDISTTOMERGEVERTEX)
		{
			f->pVF = fvv.pVF;	f->pV = fvv.pV;	break;
		}
	}
}

/* 判断下面一种情况：如果有效相交段数为0，相交点只有一个，而且这个点是角点，那么可以判断这个文化信息和这块地景不相交。*/
bool DetectIntersectStatusOne()
{
	/* 如果有效相交段个数不等于0，直接返回false。*/
	if(validBorderSegment.size() > 0)	return false;
	
	/* 统计相交点个数，大于1个，返回false */
	int intSum = 0;
	FindValidVertex fvv;
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
	{
		if (q->newIdx == 1) 
		{
			intSum++;	if(intSum > 1)	return false;
			fvv.pV = q->pV; fvv.pVF = q->pVF; fvv.pT = q->pT;
		}
	}
	
	/* 判断这个唯一的相交点是不是角点。 如果不是,返回false */
	int iFlag = 0;
	osg::Vec3 I = fvv.pVF;
	for(unsigned int i=0; i<4; i++)
	{
		double dist = LenOf2V(I, conner[i].pVF);
		if(dist < MINDISTTOMERGEVERTEX)	{ iFlag = 1; break;}
	}
	if(!iFlag) return false;
	
	/* 上述条件都满足，返回 true */
	return true;
}

/* 在两个点连成的线段中间，通过逐步的中值插入法，找到一个距离第一个点小于10米的插入点。*/
osg::Vec3 getIntersectVertexLessTenMeter(osg::Vec3 src, osg::Vec3 tgt)
{
	osg::Vec3 mid;	mid.x() = (src.x() + tgt.x()) / 2.0;	mid.y() = (src.y() + tgt.y()) / 2.0;	mid.z() = (src.z() + tgt.z()) / 2.0;
	double length = LenOf2V(src, mid);
	if(length < 10.0)
		return mid;
	else
		return getIntersectVertexLessTenMeter(src, mid);
}

/* 检查线环pm 和 物体gd 的相交情况，并增加相交点到线环里面 */
/* 返回的线环一定要是在海平面投影的线环. */
void addIntersectVertex(osg::ref_ptr<osg::Geode> gTerrain, osg::ref_ptr<osg::Geode> gCulture, osg::ref_ptr<osg::Vec3Array> pm, osg::ref_ptr<osg::Vec3Array> newpm, string pmName)
{
	if (pm->size()<3) return;

	/* 把pm里面各个点先保存到 totalVertices 里面 */
	for( osg::Vec3Array::iterator r = pm->begin(); r != pm->end(); r++ )
		*r = findPMByFeatureVertex(*r, 1, pmName);

	std::vector<FindValidVertex> fPerimeter;		//记录一个投影到海平面的外围环
	std::vector<FindValidVertex> fPerimeterNoFlat;	//记录一个原来的外围环
	std::vector<LineSeg> intersectLS;				//记录有相交的边
	std::vector<LineSeg> intersectLS_NoFlat;		//记录有相交的边

	intsResult.clear();
	fPerimeter.clear();
	fPerimeterNoFlat.clear();
	removeValidSegmentWhenWithMultiIntersections.clear();
	multiDc.clear();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 查看有几个角点在约束中
	int ic;
	int icc=0;
	for(ic=0;ic<4;ic++)
	{
		if(conner[ic].isValid > 0) {conner[ic].isValid = 1;	icc++; }
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* 找出0,3横边线上的两端坐标，利用arctan计算这条线与水平线的夹角，并求出旋转矩阵 */
	double alpha = 0;
	if (conner[3].pVF[1]-conner[0].pVF[1] != 0)
		alpha = atan((conner[3].pVF[0]-conner[0].pVF[0]) / (conner[3].pVF[1]-conner[0].pVF[1]));
	osg::Matrix _mr = osg::Matrix::rotate(alpha, osg::Z_AXIS);
	osg::Vec3  r_conner[4];
	for (ic=0;ic<4;ic++)
		r_conner[ic] = conner[ic].pVF * _mr;

	int numFi = 0;

	/* 使pm线环保证首尾相连接 */
	if (!pm->empty() && pm->front()!=pm->back())
		pm->push_back(pm->front());

	/* 遍历这个pm线环，查看线环上各个点跟地景面gd的相交关系，如果有相交，则表示该点能映射到gd面上。  */
	for( osg::Vec3Array::iterator p = pm->begin(); p != pm->end(); p++ )
	{
		/* 如果相交，表示这个点能垂直映射到 gTerrain 面上。	*/
		FindValidVertex fvv;

		/* 找到对应的投影点 */
		osg::Vec3 pFlat = getProjectionVertex(*p);

		osg::Vec3 isect;
		if( TerrainProcessingUtils::getIntersect( gTerrain, pFlat, isect ) != false )
		{
			fvv.isValid = 1;
		}
		else
		{
			//////////////////////////////////////////////////////////////////////////////////////
			// 处理两点之间的距离太大并跨过地形，不能正确判断相交的情况，

			if (numFi>0 && fPerimeter[numFi-1].isValid==0)
			{
				osg::Vec3 A, B;

				/* 找到当前点和它相邻的之前的一个点  */
				A = pFlat;
				B = fPerimeter[numFi-1].pV;

				/* 通过矩阵旋转，将空间两个点转换为平面直角坐标系上的两个点，以便检查这两个点在平面直角坐标系上面的象限情况  */
				osg::Vec3 rA   = A * _mr;
				osg::Vec3 rB   = B * _mr;

				/* 检查两个点是否在同一象限中, 如果是同一象限，则不处理 */
				bool bProcess = false;
				for (ic=0;ic<4;ic++)
				{
					/* 原点移动到 rorg 点，即以四个角点分别作为原点 */
					osg::Vec3 roA  = rA - r_conner[ic];
					osg::Vec3 roB  = rB - r_conner[ic];

					/* 如果直角坐标系内两个点的x坐标符号不同 或者y坐标符号不同，则这两个点肯定不在同一个象限内  */
					if ( (roA[0]>0 && roB[0]<0) || (roA[0]<0 && roB[0]>0) || (roA[1]>0 && roB[1]<0) || (roA[1]<0 && roB[1]>0) )
						bProcess = true;

					if (bProcess) break;
				}

				/* 如果两点不在同一象限中, 则在A,B点之间进行插值 */
				if (bProcess)
				{
					/* 插值，平均插入10个点  */
					std::vector<osg::Vec3> vec_ab;
					interpolation(A, B, vec_ab, INSERT_SEGMENTS);

					/* 求这10个点与地景面 gTerrain 的相交情况，nB表示起始相交点的序号，nE表示最终相交点的序号  */
					int nB = -1, nE = -1;
					int i = 0;
					for( std::vector<osg::Vec3>::iterator pp= vec_ab.begin(); pp != vec_ab.end(); pp++ )
					{
						if( TerrainProcessingUtils::getIntersect( gTerrain, *pp, isect ))
						{
							if (nB<0) nB = i;
							else nE = i;
						}
						else if (nE>=0)
						{
							break;
						}
						i++;
					}

					/* 如果插值序列vec_ab与地景面gd有相交点  */
					if (nB>=0 || nE>=0)
					{
						/* 求出所有相交点序列的中点  */
						osg::Vec3 iVec;
						if (nB>=0 && nE>=0)
							iVec = (vec_ab[nB]+vec_ab[nE])/2;
						else if (nB>=0)
							iVec = vec_ab[nB];

						/* 把这个中点插入到外围环fPerimeter中  */
						FindValidVertex fvv_n, fvv_nF;
						fvv_nF.pV = iVec;  fvv_nF.idx = numFi; fvv_nF.isValid = 1;
						fPerimeter.push_back(fvv_nF);

						/* 求出由原有值计算出的插值点 */
						A = *p;
						B = fPerimeterNoFlat[numFi-1].pV;
						/* 插值，平均插入10个点  */
						std::vector<osg::Vec3> vec_ab;
						interpolation(A, B, vec_ab, INSERT_SEGMENTS);
						if (nB>=0 && nE>=0)
						{
							if(nE > (int)(vec_ab.size() - 1)) nE = vec_ab.size() - 1;
							iVec = (vec_ab[nB]+vec_ab[nE])/2;
						}
						else if (nB>=0)
							iVec = vec_ab[nB];
						/* 把这个中点插入到外围环fPerimeterNoFlat中  */
						fvv_n.pV = iVec;  fvv_n.idx = numFi; fvv_n.isValid = 1; numFi++;
						fPerimeterNoFlat.push_back(fvv_n);

						/* 把顶点对加入到totalVertices里面 */
						addVertexPairIntoTotal(fvv_nF.pV, fvv_n.pV);
					}
				}
			}
			//////////////////////////////////////////////////////////////////////////////////////

			fvv.isValid = 0;
		}
		
		/* 将当前点插入到外围环fPerimeter中 */
		fvv.pV = pFlat;  fvv.idx = numFi;
		fPerimeter.push_back(fvv);

		fvv.pV = *p;  fvv.idx = numFi; numFi++;
		fPerimeterNoFlat.push_back(fvv);
	}
	
	/* 检查完毕之后，找出相交线段，isValid从0到1和从1到0变换的地方就是相交线段 */
	int i, j;
	intersectLS.clear();
	for(i=0; i< (int)(fPerimeter.size()-1); i++)
	{
		if( ((fPerimeter[i].isValid == 0) && (fPerimeter[i+1].isValid == 1)) ||
			((fPerimeter[i].isValid == 1) && (fPerimeter[i+1].isValid == 0)) ) 
		{
			LineSeg ls;
			ls.pS = fPerimeter[i].pV;   ls.idxS = fPerimeter[i].idx;	ls.validS = fPerimeter[i].isValid;
			ls.pE = fPerimeter[i+1].pV; ls.idxE = fPerimeter[i+1].idx;	ls.validE = fPerimeter[i+1].isValid;	ls.border = -1;
			intersectLS.push_back(ls);
		}
	}

	/* 找出平面点对应三维空间点的相交线段 */
	intersectLS_NoFlat.clear();
	for(i=0; i< (int)(fPerimeterNoFlat.size()-1); i++)
	{
		if( ((fPerimeterNoFlat[i].isValid == 0) && (fPerimeterNoFlat[i+1].isValid == 1)) ||
			((fPerimeterNoFlat[i].isValid == 1) && (fPerimeterNoFlat[i+1].isValid == 0)) ) 
		{
			LineSeg ls;
			ls.pS = fPerimeterNoFlat[i].pV;   ls.idxS = fPerimeterNoFlat[i].idx;	ls.validS = fPerimeterNoFlat[i].isValid;
			ls.pE = fPerimeterNoFlat[i+1].pV; ls.idxE = fPerimeterNoFlat[i+1].idx;	ls.validE = fPerimeterNoFlat[i+1].isValid;	ls.border = -1;
			intersectLS_NoFlat.push_back(ls);
		}
	}


	/* 如果有角在约束区，并且没有相交线段，并且地形块的内部点都不与文化信息相交，那么只有一种情况，这个文化信息与地形块只有一个交点，就是那个角点。*/
	if((icc > 0) && (intersectLS.size() == 0) && (AllInnerVertexNoIntersection == true)) 
	{
		/* 找到那个相交的角点，把它的onCulture状态改回来。使这个点不被归类到 invalidVertices 集合里面。*/
		for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
		{
			if(p->onCulture == 1)	p->onCulture = 0;
		}
		
		/* 退出 */
		pm->clear();
		newpm->clear();
		return;
	}

	/* 如果四个角都在约束区，并且没有相交线段，并且地形块的所有的内部点都与文化信息相交，那么就是说这块地景全部在文化信息内部*/
	/* 利用这个地块的四个角组成一个LINE_LOOP输出 */
	/* 如果这块文化信息是湖面（洋面），而且地景面的所有内部点都和这块文化信息相交，那么就表示 整个地景面都是湖面（洋面）----- 不够严肃 */
	if(((icc == 4) && (intersectLS.size() == 0) && (AllInnerVertexIntersection == true)) /* || ((pmName == "Lake") && (AllInnerVertexIntersection == true)) */ )
	{
		pm->clear();
		newpm->clear();
		for(ic=0;ic<4;ic++)
		{
			pm->push_back(connerEx[ic].pVF);
			newpm->push_back(connerEx[ic].pVF);
		}
		pm->push_back(connerEx[0].pVF);
		newpm->push_back(connerEx[0].pVF);
		
		/* 调整 checkBorderVertices 内所有的相交值。表示边上所有的点都与文化信息相交。*/
		if(	(pmName.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING) || 
				(pmName.substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) )
		{
			/* 重新确定输出线环 */
			pm->clear();
			newpm->clear();
			for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
			{
				pm->push_back(p->pVF);
				newpm->push_back(p->pVF);
			}
			pm->push_back(terrainBorderVertices[0].pVF);
			newpm->push_back(terrainBorderVertices[0].pVF);
			
			/* */			
			for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
			{
				p->onCulture = 1;
			}
			needDoSkirt = false;
		}
		subIntersectCt = faceVertice.size();
		multiDc.push_back(newpm);
		return;
	}

	/* 如果所有的边界段都和文化信息没有交点，那就直接退出。这里少考虑了一种情况：如果文化信息面完全在地景面内。*/
	if((borderHalfIntersectLS.size() == 0) && (intersectLS.size() == 0))
	{
		if(fPerimeter[0].isValid == 0)
		{			
			pm->clear();
			newpm->clear();
		}else{
			newpm->clear();
			for( osg::Vec3Array::iterator r = pm->begin(); r != pm->end(); r++ )
			{
 				*r = findPMByFeatureVertex(*r, 0, pmName);
				newpm->push_back(*r);
			}
			multiDc.push_back(newpm);
		}
		return;
	}

	/* 检查这些相交线段与边线的相交情况，结果保存在 intsResult 里面 */
	validBorderSegment.clear();
	mergeVertex.clear();
	for(j=0; j<(int)(intersectLS.size()); j++)
	{
		// 用dif保证一条线段跟四条边一定有一个交点		
		bool bInter = false;
		int dif = 0;
		while (bInter == false)
		{
			for(i=0; i<4; i++)			//i = intersectLS[j].border;		
			{
				osg::Vec3 I;

#if 1
				/* 要使 I 落在intersectLS线上，因为这是文化信息的线。*/
				intersect(intersectLS[j].pS, intersectLS[j].pE, border[i].pS, border[i].pE, I);
#else
				/* 要使 I 落在border线上，因为这是地形边界的线。*/
				intersect(border[i].pS, border[i].pE, intersectLS[j].pS, intersectLS[j].pE, I);
#endif

				//判断I是否在线段intersectLS[j]上
				if(   (((I[0] >= intersectLS[j].pS[0]-dif)&&(I[0] <= intersectLS[j].pE[0]+dif)) || ((I[0] <= intersectLS[j].pS[0]+dif)&&(I[0] >= intersectLS[j].pE[0]-dif)))
					&&(((I[1] >= intersectLS[j].pS[1]-dif)&&(I[1] <= intersectLS[j].pE[1]+dif)) || ((I[1] <= intersectLS[j].pS[1]+dif)&&(I[1] >= intersectLS[j].pE[1]-dif)))
					&&(((I[0] >= border[i].pS[0]-dif)&&(I[0] <= border[i].pE[0]+dif)) || ((I[0] <= border[i].pS[0]+dif)&&(I[0] >= border[i].pE[0]-dif)))
					&&(((I[1] >= border[i].pS[1]-dif)&&(I[1] <= border[i].pE[1]+dif)) || ((I[1] <= border[i].pS[1]+dif)&&(I[1] >= border[i].pE[1]-dif))) ) 
				{
					/* 如果交点I 坐落在相交线段intersectLS之间，而不是两个端点中的任何一个，那么就保存这个交点。*/
					FindValidVertex ir;

					/* 先判断 I 和 intersectLS[j] 两个端点间的距离，如果距离太近，可以合并成一个 */
					double distIS = LenOf2V(I, intersectLS[j].pS);
					double distIE = LenOf2V(I, intersectLS[j].pE);
					if(dif > 0)
					{
						if(distIS < 10.0)
						{
							I = intersectLS[j].pS;
						}
						if(distIE < 10.0)
						{
							I = intersectLS[j].pE;
						}
					}
					
					/* 如果交点I 坐落在相交线段intersectLS之间，而不是两个端点中的任何一个，那么就保存这个交点。*/
					if((I != intersectLS[j].pS) && (I != intersectLS[j].pE))
					{
						ir.pV = I; ir.pVF = I; ir.idx = intersectLS[j].idxS; ir.newIdx = intersectLS[j].idxE; ir.isValid = i; // 这里ir.isValid表示相交点在哪条边线上
						intsResult.push_back(ir);
					}else{
						/* 如果交点 I 和相交线段intersectLS的两个端点之一重合，那么要修改这个端点的状态。*/
						if(I == intersectLS[j].pS)
						{
							if(intersectLS[j].validS == 0)
							{
								intersectLS[j].validS = 1;
								fPerimeter[intersectLS[j].idxS].isValid = 1;
							}
						}
						if(I == intersectLS[j].pE)
						{
							if(intersectLS[j].validE == 0)
							{
								intersectLS[j].validE = 1;
								fPerimeter[intersectLS[j].idxE].isValid = 1;
							}
						}
					}
					
					/* 根据这个投影点，计算出与其对应的空间点。*/
					/* 首先转换成经纬高坐标系 */
					osg::Vec3 I2;
					int idxPv;
					if((idxPv = FindVertexFromIntersect(I)) == -1)
					{
						Geodesy geoI, geoS, geoE;
						Carton  cat;
						cat.x = I.x();	cat.y = I.y(); cat.z = I.z();	CartToGeod(&cat, &geoI);
						osg::Vec3 S; S = intersectLS_NoFlat[j].pS;
						cat.x = S.x();	cat.y = S.y(); cat.z = S.z();	CartToGeod(&cat, &geoS);
						osg::Vec3 E; E = intersectLS_NoFlat[j].pE;
						cat.x = E.x();	cat.y = E.y(); cat.z = E.z();	CartToGeod(&cat, &geoE);
						/* 根据 空间直线两点式 */
						if (geoE._elevation - geoS._elevation == 0)
						{
							geoI._elevation = geoE._elevation;
						}
						else if (geoE._lon - geoS._lon != 0)
						{
							//(x-x1)/(x2-x1)=(z-z1)/(z2-z1)
							geoI._elevation = (geoE._elevation - geoS._elevation)*(geoI._lon-geoS._lon)/(geoE._lon - geoS._lon) + geoS._elevation;
						}
						else
						{
							// (y-y1)/(y2-y1)=(z-z1)/(z2-z1)
							geoI._elevation = (geoE._elevation - geoS._elevation)*(geoI._lat-geoS._lat)/(geoE._lat - geoS._lat) + geoS._elevation;
						}
						GeodToCart(&geoI, &cat);	I2.set(cat.x, cat.y, cat.z);
						addVertexIntoIntersect(I2);
					}else
					{
						I2 = IntersectVertices[idxPv].pV;
					}
					addVertexIntoTotal(I2);

					/* 把相交点加入到有效边上点的集合 */
					ir.pV = I; ir.idx = 0; ir.newIdx = 0; ir.isValid = 0; // 这里ir.isValid表示相交点在哪条边线上
					validBorderVertices.push_back(ir);

					ir.pVF = I; ir.pV = I2; ir.idx = 0; ir.newIdx = 0; ir.isValid = 0; // 这里ir.isValid表示相交点在哪条边线上
					switch(i)
					{
					case 0:
						ir.pT.set(0.0, 0.555555555);
						break;
					case 1:
						ir.pT.set(0.555555555, 1.0);
						break;
					case 2:
						ir.pT.set(1.0, 0.555555555);
						break;
					case 3:
						ir.pT.set(0.555555555, 0.0);
						break;
					}
					allBorderVertices.push_back(ir);
					
					/* 找到这个点在 哪个 borderIntersectLS 里面 */
					unsigned int k, index = -1;
					double minDistant = 999999999.999;
					for(k=0; k<borderIntersectLS.size(); k++)
					{
						if(borderIntersectLS[k].border != i) continue;
						double dist = LenOf2V(I, borderIntersectLS[k].pS) + LenOf2V(I, borderIntersectLS[k].pE);
						if(	dist < minDistant)
						{
							minDistant = dist;
							index = k;
						}
					}
					
					if(index != -1)
					{
						int iFlag = -1;			//这个标志标识ls.pE 是取的是边界段的 pS 点 还是 pE 点。
						int addValidBorderFlag = 1;
						LineSeg ls;
						ls.pS = I;   ls.idxS = -1;	ls.validS = 1;
						if(( borderIntersectLS[index].validS == 1) && ( borderIntersectLS[index].validE == 0))
						{	iFlag  = 0; ls.pE = borderIntersectLS[index].pS; ls.idxE = borderIntersectLS[index].idxS;	ls.validE = 1;	}	//ls.validE 判断这个边界段的两个端点与文化信息相交情况，2:都相交；1:只有一个相交；0:都不相交
						if(( borderIntersectLS[index].validS == 1) && ( borderIntersectLS[index].validE == 1))
						{	iFlag  = 0; ls.pE = borderIntersectLS[index].pS; ls.idxE = borderIntersectLS[index].idxS;	ls.validE = 2;	}
						if(( borderIntersectLS[index].validS == 0) && ( borderIntersectLS[index].validE == 1))
						{	iFlag  = 1; ls.pE = borderIntersectLS[index].pE; ls.idxE = borderIntersectLS[index].idxE;	ls.validE = 1;	}
						if(( borderIntersectLS[index].validS == 0) && ( borderIntersectLS[index].validE == 0))
						{	iFlag  = 1; ls.pE = borderIntersectLS[index].pE; ls.idxE = borderIntersectLS[index].idxE;	ls.validE = 0;	}
						ls.border = index;		//border表示当前相交段是第几段

						/* 从半相交边界段集合中，移除与本次相交结果相同的段 */
						if(ls.validE == 1)
						{
							for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
							{
								if((r->idxS == ls.idxE) && (iFlag == 0)) r->border = -2;
								if((r->idxE == ls.idxE) && (iFlag == 1)) r->border = -2;
							}
						}
						
						/* 计算这个交点和边界段上的两个端点之间的距离，如果距离<1.0，可认为该交点和端点是一个点*/
						/* 注意：这里只判断xy,不判断z，为了和三角化算法保持一致。Galen:2013-02-06 */
						double distS = LenOf2VF(I, borderIntersectLS[index].pS);
						double distE = LenOf2VF(I, borderIntersectLS[index].pE);
						
						/* 如果这个交点所在边的两个端点，都是和文化信息面相交的点 （即ls.validE = 2），并且这个交点不和任意一个端点重合， */
						/* 那就要判断这个交点的两边，到底是哪一边和文化信息有相交？ */
						/* 目前简单地采用插入点相交判断法。*/
						if((ls.validE == 2) && (distS >= MINDISTTOMERGEVERTEX) && (distE >= MINDISTTOMERGEVERTEX))
						{
							/* 计算交点与后一个端点之间的插入点是否与文化信息相交。*/
							osg::Vec3 IntVer = getIntersectVertexLessTenMeter(I, borderIntersectLS[index].pE);
							osg::Vec3 isect;
							if( TerrainProcessingUtils::getIntersect( gCulture, IntVer, isect ) != false )
							{
								ls.pE = borderIntersectLS[index].pE;	ls.idxE = borderIntersectLS[index].idxE;
							}
							removeValidSegmentWhenWithMultiIntersections.push_back(index);
						}
						
						/* 如果当前交点与S边的端点接近一致。并且当前边界段两端的状态相同。 */
						if(distS < MINDISTTOMERGEVERTEX) 
						{
							/* 这个变量表示不加入有效相交段 validBorderSegment */
							addValidBorderFlag = 0;

							/* 替换faceVertice 中的相关点。*/
							ModifyFaceVertexWhenMerge(ir);
							
							/* 找到 checkBorderVertices 中相同的点并处理 */
							FindVertexInCBVsets(I);
							
							borderIntersectLS[index].pS = I;						//既然可以当作一个点，那么就把边界段的端点换成这个相交点。
							if (borderIntersectLS[index].validS == 0)	borderIntersectLS[index].validS = 1;
							int idx = index - 1; if(idx <0) idx += borderIntersectLS.size();
							borderIntersectLS[idx].pE = I;				//既然可以当作一个点，那么就把边界段的端点换成这个相交点。
							if (borderIntersectLS[idx].validE == 0)	borderIntersectLS[idx].validE = 1;
							
							/* 判断边界顶点集terrainBorderVertices 里面这个位置的原来相交判断值newIdx是否为1 */
							int idxTBV = borderIntersectLS[index].idxS;		//这个值是terrainBorderVertices里面对应下标。
							if(terrainBorderVertices[idxTBV].newIdx >= 0)
							{
								/* 判断当前index的borderIntersectLS边界相交段，和index-1的边界相交段，如果这两个相交段都是没有与文化信息相交相交，*/
								/* 说明当前文化信息穿过 index 相交段内但是与边界点无相交。这个时候，要加入validBorderSegment相交段，以便后面的计算。*/
								/* idx = (index - 1)%borderIntersectLS.size() */
								if((borderIntersectLS[idx].validS == 0) && (borderIntersectLS[index].validS == 0) && (borderIntersectLS[index].validE == 0))
								{
									addValidBorderFlag = 1;
								}

								/* 判断这个下标的前一个值是否和文化信息相交（newIdx == 1），*/
								/* 如果前面那个值也是相交，那么这一个线段原本是个半相交的线段，更改之后名义上成了全相交的线段，但它还是半相交的线段，所以要预先排出它。*/
								if((idxTBV != 0) && (terrainBorderVertices[idxTBV].newIdx == 0))
								{
									int idxPrev = idxTBV - 1;
									if(terrainBorderVertices[idxPrev].newIdx == 1)
										removeValidSegmentWhenWithMultiIntersections.push_back(idxPrev);
								}

								/* 如果这个下标值是0，那要把数组最后面那个值也赋成一样的值，因为这是一个线环。*/
								if(idxTBV == 0)
								{
									int ii = terrainBorderVertices.size() - 1;

									int idxPrev = ii - 1;
									if((terrainBorderVertices[idxPrev].newIdx == 1) && (terrainBorderVertices[ii].newIdx == 0))
										removeValidSegmentWhenWithMultiIntersections.push_back(idxPrev);

									/* 用新值代替旧值，且标记这个点肯定和文化信息相交（因为这个点也是文化信息边界的点）。*/
									terrainBorderVertices[ii].pVF = I;
									terrainBorderVertices[ii].newIdx = 1;
								}

								/* 用新值代替旧值，且标记这个点肯定和文化信息相交（因为这个点也是文化信息边界的点）。*/
								terrainBorderVertices[idxTBV].pVF = I;
								terrainBorderVertices[idxTBV].newIdx = 1;		// 1代表这个点同文化信息相交，0：代表这个点同文化信息不相交。
							}
							
							/* 相交点是角点的处理 */
							CheckSpecialIntVertex(I);
						
							/* 还要查找先前做好的validBorderSegment里面有没有和I重合的点，如果有，就替换过来 */
							for(unsigned int m = 0; m<validBorderSegment.size(); m++)
							{
								double distVS = LenOf2V(I, validBorderSegment[m].pS);
								double distVE = LenOf2V(I, validBorderSegment[m].pE);
								if(distVS < MINDISTTOMERGEVERTEX) validBorderSegment[m].pS = I;
								if(distVE < MINDISTTOMERGEVERTEX) validBorderSegment[m].pE = I;
							}
							
							/* 保存合并点的序号 */
							int iFlag = 0;
							for(unsigned int k = 0;k<mergeVertex.size();k++)
							{
								if(index == mergeVertex[k]) {	iFlag = 1; break; }
							}
							if(!iFlag)	mergeVertex.push_back(index);
						}

						/* 如果当前交点与E边的端点接近一致。并且当前边界段两端的状态相同。 */
						if(distE < MINDISTTOMERGEVERTEX) 
						{
							/* 替换faceVertice 中的相关点。*/
							ModifyFaceVertexWhenMerge(ir);
							
							/* 找到 checkBorderVertices 中相同的点并处理 */
							FindVertexInCBVsets(I);

							borderIntersectLS[index].pE = I;						//既然可以当作一个点，那么就把边界段的端点换成这个相交点。
							if (borderIntersectLS[index].validE == 0)	borderIntersectLS[index].validE = 1;
							unsigned int idx = index + 1; if(idx >= borderIntersectLS.size()) idx = 1;
							borderIntersectLS[idx].pS = I;				//既然可以当作一个点，那么就把边界段的端点换成这个相交点。
							if (borderIntersectLS[idx].validS == 0)	borderIntersectLS[idx].validS = 1;

							/* 判断边界顶点集terrainBorderVertices 里面这个位置的原来相交判断值newIdx是否为1 */
							int idxTBV = borderIntersectLS[index].idxE;		//这个值是terrainBorderVertices里面对应下标。
							if(terrainBorderVertices[idxTBV].newIdx >= 0)
							{
								/* 用新值代替旧值，且标记这个点肯定和文化信息相交（因为这个点也是文化信息边界的点）。*/
								terrainBorderVertices[idxTBV].pVF = I;
								terrainBorderVertices[idxTBV].newIdx = 1;

								/* 这里还要考虑首尾边界的问题。*/
								/* 因为terrainBorderVertices最后一个点和第一个点是一个点。*/
								if(idxTBV == terrainBorderVertices.size() -1 )
								{
									terrainBorderVertices[0].pVF = I;
									terrainBorderVertices[0].newIdx = 1;
								}
							}
							addValidBorderFlag = 0;

							/* 相交点是角点的处理 */
							CheckSpecialIntVertex(I);
						
							/* 还要查找先前做好的validBorderSegment里面有没有和I重合的点，如果有，就替换过来 */
							for(unsigned int m = 0; m<validBorderSegment.size(); m++)
							{
								double distVS = LenOf2V(I, validBorderSegment[m].pS);
								double distVE = LenOf2V(I, validBorderSegment[m].pE);
								if(distVS < MINDISTTOMERGEVERTEX) validBorderSegment[m].pS = I;
								if(distVE < MINDISTTOMERGEVERTEX) validBorderSegment[m].pE = I;
							}

							/* 保存合并点的序号 */
							int iFlag = 0;
							for(unsigned int k = 0;k<mergeVertex.size();k++)
							{
								if((index+1) == mergeVertex[k]) {	iFlag = 1; break; }
							}
							if(!iFlag)	mergeVertex.push_back(index+1);
						}
						
						if(addValidBorderFlag)	validBorderSegment.push_back(ls);
					}

					/* 设置标记 */
					bInter = true;
				}
			}
			dif++;
		}
	}//for(j)
	
	/* 	判断下面一种情况：如果有效相交段数为0，相交点只有一个，而且这个点是角点，
		那么可以判断这个文化信息和这块地景不相交。*/
	if(DetectIntersectStatusOne())
	{
		multiDc.clear();
		return;
	}

	/* 检测特殊情况 */
	CheckSpecialStatus();

	/* 这里要处理一种情况，如果文化信息曲线与一段边界有多次的相交，该怎么办？*/
	CheckSegmentsIntersection();

	/* 把边线上其他的有效边也加入到validBorderSegment里面 */
	std::vector<LineSeg> tmpSegment;
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
		tmpSegment.push_back(*q);

	int borderIdx = 0;
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
	{
		if(q != terrainBorderVertices.begin())
		{
			std::vector<FindValidVertex>::iterator p = q - 1;
			if( (p->newIdx == 1) && (q->newIdx == 1) )
			{
				LineSeg ls;
				ls.pS = p->pVF;   ls.idxS = p->idx;	ls.validS = p->newIdx;
				ls.pE = q->pVF;   ls.idxE = q->idx;	ls.validE = q->newIdx;	ls.border = p->idx;

				/* 检查这对线段之间，是否已经存在一个相交边界线段*/
				/* 检查方法是：判断 validBorderSegment 里面的 border，如果当前段和border段一致，说明已经存在一条相交边界段了。*/
				int jFlag = 0;
				for( std::vector<LineSeg>::iterator r= tmpSegment.begin(); r != tmpSegment.end(); r++ )
				{
					if(r->border == borderIdx) { jFlag = 1; break;	}
				}
				if(jFlag) goto nextBorder;
				
				/* 检查这个有效段，是否中间有多个交点，如果有的话，则不增加这个有效段。 */
				int iFlag = 0;
				for(unsigned int i=0; i<removeValidSegmentWhenWithMultiIntersections.size(); i++)
				{
					if(ls.idxS == removeValidSegmentWhenWithMultiIntersections[i]) { iFlag = 1; break; }
				}
				if(!iFlag)	validBorderSegment.push_back(ls);
			}
nextBorder:
			/* 指向下一段 */
			borderIdx++;
		}
	}
	
	/* 把这些相交点插进pm线环里头 */
	/* Step1:先把intsResult 按照idx序号从小到大顺序排列 */
	for(i=0; i<(int)(intsResult.size()-1); i++)
	{
		for(j=i+1; j<(int)(intsResult.size()); j++)
		{
			if(intsResult[i].idx > intsResult[j].idx)
			{
				FindValidVertex tmp;
				tmp.pV = intsResult[i].pV; tmp.idx = intsResult[i].idx; tmp.newIdx = intsResult[i].newIdx;
				intsResult[i].pV = intsResult[j].pV; intsResult[i].idx = intsResult[j].idx; intsResult[i].newIdx = intsResult[j].newIdx;
				intsResult[j].pV = tmp.pV; intsResult[j].idx = tmp.idx; intsResult[j].newIdx = tmp.newIdx;
			}
		}
	}

	/* Step2:把排序后的intsResult，插入到fPerimeter序列里面  */
	for(i=intsResult.size()-1; i>=0; i--)
	{
		j = fPerimeter.size();
		FindValidVertex tmp;
		fPerimeter.push_back(tmp);
		while(j > intsResult[i].newIdx) { fPerimeter[j].pV = fPerimeter[j-1].pV;  fPerimeter[j].isValid = fPerimeter[j-1].isValid; fPerimeter[j].idx = fPerimeter[j-1].idx; fPerimeter[j].newIdx = j; j--;}
		fPerimeter[intsResult[i].newIdx].pV = intsResult[i].pV;
		fPerimeter[intsResult[i].newIdx].isValid = intsResult[i].isValid + 2;  // 用2来表示该点是新插入的相交点，并表示相交点在哪条边上，2表示边0，3表示1边...
		fPerimeter[intsResult[i].newIdx].idx = -1;fPerimeter[intsResult[i].newIdx].newIdx = intsResult[i].newIdx;
	}

	/* 把fPerimeter线环里面的有效线段也都添加到 validBorderSegment 里面来 */
	for( std::vector<FindValidVertex>::iterator q= fPerimeter.begin(); q != fPerimeter.end(); q++ )
	{
		if(q != fPerimeter.begin())
		{
			std::vector<FindValidVertex>::iterator p = q - 1;
			if( (p->isValid > 0) && (q->isValid > 0) )
			{
				LineSeg ls;
				ls.pS = p->pV;   ls.idxS = p->idx;	ls.validS = p->newIdx;
				ls.pE = q->pV;   ls.idxE = q->idx;	ls.validE = q->newIdx;	ls.border = -1;
				validBorderSegment.push_back(ls);
			}
		}
	}

	/* 如果没有有效边段，那么直接退出去。*/
	if(validBorderSegment.size() == 0)
	{
		pm->clear();
		newpm->clear();
		return;
	}


	/* 把 validBorderSegment 组成一个线环 */
	osg::ref_ptr<osg::Vec3Array> pmLL = new osg::Vec3Array;
	CreateLineLoopFromValidBorderSegment(pmLL);
	
	/* 如果 validBorderSegment 没有形成有效线环。*/
	if(multiDc.size() == 0)	return;
	if(multiDc[0]->getNumElements() <= 2)
	{
		multiDc.clear();
		return;
	}

	/* Add by Galen 2013-02-18 */
	/* 如果当前线环的文化信息是湖，那么把这个线环上各个点先映射一遍。保存到 totalVertices 数组里面*/
	if(	(pmName.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING) || 
			(pmName.substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) )
	{
		for(std::vector< osg::ref_ptr<osg::Vec3Array> >::iterator mq = multiDc.begin(); mq != multiDc.end(); mq++)
		{
			osg::ref_ptr<osg::Vec3Array> currVa = *mq;
			for(osg::Vec3Array::iterator q = currVa->begin(); q != currVa->end(); ++q)
			{
				addVertexIntoTotalWithElevation(*q, elevationOfLake);
			}
		}
	}
	/* Add by Galen End */


	/* 将找到的相交点，全部加入到有效点的判断集合里面。 */
	for(i=0; i<(int)(intsResult.size()); i++)
	{
		intersectlist.push_back(intsResult[i]);
	}

	for (ic=0;ic<4;ic++)
	{
		if (conner[ic].isValid == 1) 
			intersectlist.push_back(conner[ic]);
	}
	
}


/* 检查va,vb,vc三点所组成的三角形，是否与地形方块terrain相交？
   可以重新写一个算法：
   首先计算地景面的包围球，
   再计算每个三角面的包围球，
   只要两个球心之间的距离大于两球半径之和，就说明面不相交。返回0。
   如果两个球心之间的距离小于两球半径之和，那么再用插值的办法来仔细判断。
   做完上述判断后，再用德罗尼算法的 handleOverlaps()函数，求出地景面LINE_LOOP和三角形LINE_LOOP的公共交点个数。
   就不需要再做插值算法了。
 */
int checkTriwithTerrain(osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vc, osg::Geode *terrain, osg::ref_ptr<osg::Vec3Array> tbLoop)
{
	const osg::BoundingSphere bsTerrain = terrain->getBound();
	osg::ref_ptr<osg::Geode> tri = new osg::Geode;

	/* 首先，包围球过滤 */	
	/* 把这三个点定义成一个LINE_LOOP的Geode，再计算出包围球 */
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
	osg::ref_ptr<osg::Vec3Array> v = new osg::Vec3Array;
	v->push_back(getProjectionVertex(va));
	v->push_back(getProjectionVertex(vb));
	v->push_back(getProjectionVertex(vc));
	v->push_back(getProjectionVertex(va));
	geom->setVertexArray(v);
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, v->size()));
	tri->addDrawable(geom);
	const osg::BoundingSphere bsTriangle = tri->getBound();

	/* 计算两个包围球球心之间的距离，和两个包围球半径之和 */	
	float sumRadius = bsTerrain.radius() + bsTriangle.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsTriangle.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return 0;
	

	/* 其次，顶点过滤 */
	/* 三角形三个顶点中，如果有任意一个顶点同地形面相交，则返回1。*/
	std::vector<osg::Vec3> vec_abc;
	vec_abc.push_back(va);
	vec_abc.push_back(vb);
	vec_abc.push_back(vc);
	for( std::vector<osg::Vec3>::iterator p= vec_abc.begin(); p != vec_abc.end(); p++ )
	{
		osg::Vec3 isect;
		if( TerrainProcessingUtils::getIntersect( terrain, *p, isect ) != false )
			return 1;
	}

	/* 再次，德罗尼求交过滤 */	
	/* 如果三角形的三个顶点与地景面都没有交点，那么有可能三条边与地景面有交点。以下运算就是计算三条边与地景面交点的个数。*/
	osg::ref_ptr<ArealConstraint> dcTriangle = new ArealConstraint;
	dcTriangle->setVertexArray( v.get() );
	dcTriangle->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, v->size()) );
	unsigned int sum1 = v->size();
	osg::ref_ptr<ArealConstraint> dcTerrain = new ArealConstraint;
	dcTerrain->setVertexArray( tbLoop.get() );
	dcTerrain->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, tbLoop->size()) );
	unsigned int sum2 = tbLoop->size();
	osg::ref_ptr<ArealConstraint> dcMerge = new ArealConstraint();
	dcMerge->merge(dcTerrain.get());
	dcMerge->merge(dcTriangle.get());

	dcMerge->handleOverlaps1();
	unsigned int sum = (dcMerge->getVertexArray())->getNumElements();

	if(sum > (sum1 + sum2 + 2))
		return 1;
	else
		return 0;
}

/* 检查裙边的三个点是否都是有效点，只要有一个不是有效点，就返回0。*/
int checkTriOnSkirt(osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vc)
{     
	/* 检查这三个点是否都是有效点，只要有一个不是有效点，就返回0。*/
	for( std::vector<osg::Vec3>::iterator p= invalidVertices.begin(); p != invalidVertices.end(); p++ )
	{
		if (va==*p || vb==*p || vc==*p)
			return 0;
	}
	return 1;
}


/* 检查面上的三个点。*/
/* 如果三个点中, 有两个点在文化信息上, 而且这两个点不邻接, 那么这个三角面是非法的。 */
struct vIdx
{
	osg::Vec3 Vtx;
	int idx;
};
int checkTriOnDclist(osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vc)
{
	return 1;
	
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		int index = 0;
		osgUtil::DelaunayConstraint * dc = dynamic_cast<osgUtil::DelaunayConstraint *> (p->get());
		if(dc->getNumPrimitiveSets() > 0)
		{
			std::vector<vIdx> allvtxs;	allvtxs.clear();
			osg::ref_ptr<osg::Vec3Array> vas = dynamic_cast<osg::Vec3Array *> (dc->getVertexArray());
			for( osg::Vec3Array::iterator r = vas->begin(); r != vas->end(); r++ )
			{
				osg::Vec3 v = *r;
				if( v == va )	{	vIdx vi;	vi.Vtx = va;	vi.idx = index; allvtxs.push_back(vi);	}
				if( v == vb )	{	vIdx vi;	vi.Vtx = vb;	vi.idx = index; allvtxs.push_back(vi);	}
				if( v == vc )	{	vIdx vi;	vi.Vtx = vc;	vi.idx = index; allvtxs.push_back(vi);	}
				index++;
			}

			if(allvtxs.size() == 2)
			{
				if((abs(allvtxs[0].idx - allvtxs[1].idx) != 1) && (abs(allvtxs[0].idx - allvtxs[1].idx) != (vas->size() - 1) ))
				{
					unsigned int i = 0;
					return 0;
				}
			}
			if(allvtxs.size() == 3)
			{
				if(allvtxs[0].idx != 0)
					return 0;
			}
		}
	}
	
	return 1;
}

/* 检查面上边界的三个点是否都是有效点，只要有一个不是有效点，就返回0。*/
int checkTriOnBorder(FindValidVertex *a, FindValidVertex *b, FindValidVertex *c)
{     
	osg::Vec3 va = a->pVF, vb = b->pVF, vc = c->pVF;
	osg::Vec3 sa = a->pV, sb = b->pV, sc = c->pV;
	
	/* 检查这三个点是否都是有效点，只要有一个不是有效点，就返回0。*/
	for( std::vector<osg::Vec3>::iterator p= invalidVertices.begin(); p != invalidVertices.end(); p++ )
	{
		if (va==*p || vb==*p || vc==*p)
		{
			int i=0;
			return 0;
		}
	}

	/* 不排除三个点中有两个点相同 */
	if(( LenOf2V(va,vb) < MINDISTTOMERGEVERTEX ) || ( LenOf2V(va,vc) < MINDISTTOMERGEVERTEX ) || ( LenOf2V(vc,vb) < MINDISTTOMERGEVERTEX ))
	{
		int i=0;
		return 0;
	}

	/* 检查这三个点是否都是边界上的点，如果是边界点，那么就取得对应的纹理坐标。*/
	osg::Vec2 ta, tb, tc;
	ta.set(-1.0, -1.0);	tb.set(-1.0, -1.0);	tc.set(-1.0, -1.0);	
	int iBorders = 0;
	for( std::vector<FindValidVertex>::iterator p = allBorderVertices.begin(); p != allBorderVertices.end(); p++ )
	{
		if(va==p->pVF) {	if((ta.x()==-1.0) && (ta.y() == -1.0))	{	ta = p->pT;	iBorders++;	}	};
		if(vb==p->pVF) {	if((tb.x()==-1.0) && (tb.y() == -1.0))	{	tb = p->pT;	iBorders++;	}	};
		if(vc==p->pVF) {	if((tc.x()==-1.0) && (tc.y() == -1.0))	{	tc = p->pT;	iBorders++;	}	};
	}

	/* 面积高度检验法， */
	/* 首先求出这个三角形的面积 */
	double area = getArea(sa, sb, sc);
	/* 求出最长的一边 */
	double maxLen = distanceV(sa, sb);
	double lenBC = distanceV(sb, sc);	if(maxLen < lenBC) maxLen = lenBC;
	double lenAC = distanceV(sa, sc);	if(maxLen < lenAC) maxLen = lenAC;
	double heigh = 2 * area / maxLen;
	/* 如果这个三角形的高太小，并且三角形至少有两个点在地形的边上。那这个三角形就不要了。*/
	if( heigh < 1.0 )
	{
		if(iBorders > 1)
			return 0;
		if(( heigh < 0.5 ) && ( iBorders > 0))
			return 0;
	}

	/* 如果任意一个点不是边界上的点，则返回 1。*/
	if( (ta.x() == -1.0) || (tb.x() == -1.0) || (tc.x() == -1.0) )	return 1;

	/* 如果三个点都是边界上的点，但纹理坐标有一位相同，说明这三个点在同一侧的边上，则返回0。*/
	if( ((ta.x() == tb.x()) && (ta.x() == tc.x())) || ((ta.y() == tb.y()) && (ta.y() == tc.y())) ) 
	{
		int i=0;
		return 0;
	}

	return 1;
}

#ifdef ADD_SKIRT
/* 重新制作所有的裙边面 */
std::vector<FindValidVertex> cultSkirtVertice;
std::vector< osg::ref_ptr<osg::DrawElementsUInt> > cultSkirtDrawElements;

void makeSkirt(std::vector<osg::ref_ptr<osg::Vec3Array>> &perimeterList, int currSkirtVertexIndex)
{
	Geodesy geo;
	Carton  cat;
	
	/* 增加裙边*/
	cultSkirtVertice.clear();
	cultSkirtDrawElements.clear();
	
	for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perimeterList.begin(); q != perimeterList.end(); q++ )
	{
		std::vector<FindValidVertex> fPerimeterlist;
		for( osg::Vec3Array::iterator p = q->get()->begin(); p != q->get()->end(); p++ )
		{
			FindValidVertex fv;
			fv.pVF = *p;
			fv.pV  = *p;
			fv.newIdx = -1;
			
			/* 先判断出这个点,是face面上的点（isValid = 1），还是culture的点（isValid = 0），*/
			for( std::vector<FindValidVertex>::iterator pp= faceVertice.begin(); pp != faceVertice.end(); pp++ )
			{
				fv.isValid = 0;
				if((pp->pVF == fv.pVF) || (pp->pV == fv.pV))
				{
					fv.isValid = 1;
					break;
				}
			}

			/* 找出所有合法PerimeterVec的新序号 */
			for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
			{
				//if(pp->pV == fv.pV)
				float dx = fabs(pp->pVF.x() - fv.pVF.x());
				float dy = fabs(pp->pVF.y() - fv.pVF.y());
				float dz = fabs(pp->pVF.z() - fv.pVF.z());
				if((dx <= 1.0) && (dy <= 1.0) && (dz <= 1.0))
				{
					fv.newIdx = pp->newIdx;
					fv.pV = pp->pV;
					break;
				}

				dx = fabs(pp->pV.x() - fv.pV.x());
				dy = fabs(pp->pV.y() - fv.pV.y());
				dz = fabs(pp->pV.z() - fv.pV.z());
				if((dx <= 1.0) && (dy <= 1.0) && (dz <= 1.0))
				{
					fv.newIdx = pp->newIdx;
					fv.pVF = pp->pVF;
					break;
				}
			}

			/* fv.newIdx 表示在有效点集合里面没有找到这个线环上的点。这个很罕见，但也会出现。*/
			//assert(fv.newIdx != -1);
			if(fv.newIdx != -1)
				fPerimeterlist.push_back(fv);
		}
		if (fPerimeterlist.empty()) continue;
			
		/* 逐个点检查，如果这个点的前一个点和后一个点都是culture上的点，但是这个点确是face上的点，那么调整这个点的高程。*/
		unsigned int fsize = fPerimeterlist.size();
		for(unsigned int i=0; i<fsize; i++ )
		{
			/* 当前点 */
			FindValidVertex *fvCurr; fvCurr = &fPerimeterlist[(i+fsize)%fsize];
			/*前一个点*/
			FindValidVertex *fvPrec; fvPrec = &fPerimeterlist[(i-1+fsize)%fsize];
			/*后一个点*/
			FindValidVertex *fvSucc; fvSucc = &fPerimeterlist[(i+1+fsize)%fsize];
			if((fvPrec->isValid == 0) && (fvSucc->isValid == 0) && (fvCurr->isValid == 1))
			{
				cat.x = fvPrec->pV[0]; cat.y =  fvPrec->pV[1]; cat.z =  fvPrec->pV[2]; 
				CartToGeod(&cat, &geo);
				float elevPrec = geo._elevation;

				cat.x = fvSucc->pV[0]; cat.y =  fvSucc->pV[1]; cat.z =  fvSucc->pV[2]; 
				CartToGeod(&cat, &geo);
				float elevSucc = geo._elevation;

				cat.x = fvCurr->pV[0]; cat.y =  fvCurr->pV[1]; cat.z =  fvCurr->pV[2]; 
				CartToGeod(&cat, &geo);
				float elevCurr = geo._elevation;
				geo._elevation = (elevPrec + elevSucc) / 2.0;
				GeodToCart(&geo, &cat);
				fvCurr->pV[0] = cat.x;	fvCurr->pV[1] = cat.y;	fvCurr->pV[2] = cat.z;
				
				/* 更新 validVertices 里面的对应值 */
				validVertices[fvCurr->newIdx].pV = fvCurr->pV;
			}
		}
		
		/* 新增点集合和DrawElements */
		osg::ref_ptr<osg::DrawElementsUInt> deui = new osg::DrawElementsUInt(GL_TRIANGLE_STRIP);
		osg::ref_ptr<osg::DrawElementsUInt> revdeui = new osg::DrawElementsUInt(GL_TRIANGLE_STRIP);
		for( std::vector<FindValidVertex>::iterator p = fPerimeterlist.begin(); p != fPerimeterlist.end(); p++ )
		{
			FindValidVertex fv;
			
			cat.x = p->pV[0]; cat.y =  p->pV[1]; cat.z =  p->pV[2]; 
			CartToGeod(&cat, &geo);
			geo._elevation -= skirtLength;
			GeodToCart(&geo, &cat);
			fv.pV[0] = cat.x;	fv.pV[1] = cat.y;	fv.pV[2] = cat.z;	

			fv.newIdx = currSkirtVertexIndex;
			cultSkirtVertice.push_back(fv);
			currSkirtVertexIndex++;

			deui->push_back(p->newIdx);
			deui->push_back(fv.newIdx);

			revdeui->push_back(fv.newIdx);
			revdeui->push_back(p->newIdx);
		}
		cultSkirtDrawElements.push_back(deui);
		cultSkirtDrawElements.push_back(revdeui);

	}
}
#endif	//ADD_SKIRT

/* Add by Galen 2013-02-20 */
/* 
   计算： va, vb 两点与 vO 的距离之和。
 */
double getDistanceFromTwoVertexToVertex( osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vO)
{

	/* 计算va与vO点之间的距离。*/	
	osg::Vec3 Vdisa = getProjectionVertex(vO) - getProjectionVertex(va);
	double Disa = Vdisa.length();

	/* 计算vb与vO点之间的距离。*/	
	osg::Vec3 Vdisb = getProjectionVertex(vO) - getProjectionVertex(vb);
	double Disb = Vdisb.length();

	return Disa + Disb;	
}

/* Add by Galen 2013-02-21 
   将一个Geometry面转换为边的集合
   输入：osg::Geometry *
   输出：std::vector<Edge>
 */
void getEdgesOfGeometry(std::vector<TerrainProcessingUtils::Edge> *kept, osg::Geometry *g)
{
	std::vector<TerrainProcessingUtils::Edge>   edges;
	osg::TriangleFunctor<TerrainProcessingUtils::MakeEdgeList> mel;
	mel.setEdgeVector(&edges);
	g->accept( mel );
	
	/* 判断edges里面所有的边，如果这个边出现多次，就把这个边标注为false，表示这个边是面内部的边。2012-10-10 */
	for( std::vector<TerrainProcessingUtils::Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		std::vector<TerrainProcessingUtils::Edge>::iterator q = p;
		q++;
		while(q != edges.end() )
		{
			if( *p == *q )
			{
				p->keep = false;
				q->keep = false;
			}
			q++;
		}
	}

	/* 把edges里面多次出现的边排除掉，只留下只出现一次的边，这个边一定是面的边界上的边。2012-10-10 */
	for( std::vector<TerrainProcessingUtils::Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		if(( p->keep ) && (p->A != p->B ))
		kept->push_back( *p );
	}
}	

 
/* 下面实现这样一个算法：计算一个地形面上各个点同文化信息面的边界的距离，要求：
   1、这个地形面上的点与文化信息面不相交。
   2、这个地形面上的点距离文化信息面最近的边线的距离小于一个常数值。
 */
void getVextexListNearToCulture(osg::Geode* gdT, osg::Geometry *gmC)
/* 输入：gdT, 映射后的地景面，gmC, 映射后的文化信息面 */
/* 输出：地景面上符合条件的点集合 */
{
	/* Add by Galen 2013-02-22 */
	/* 首先判断几何体g整体的包围球是否与地形gd包围球有重合 */
	const osg::BoundingSphere bsTerrain  = gdT->getBound();
	const osg::BoundingSphere bsLakeGmtry = gmC->getBound();

	/* 计算两个包围球球心之间的距离，和两个包围球半径之和 */	
	float sumRadius = bsTerrain.radius() + bsLakeGmtry.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsLakeGmtry.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return ;
	/* Add by Galen 2013-02-22 End */
	
	
	/* 创建一个文化信息的Geode, 用来计算相交。*/
	osg::ref_ptr<osg::Geode> flatFeatureGd = new osg::Geode;
	flatFeatureGd->addDrawable(gmC);
	
	/* 将这个几何体转换为线环，然后将线环转换为边集合。*/
	/* Add by Galen 2013-02-22 先判断这个湖的几何体是否已经边段化。*/
	osg::ref_ptr<osg::Vec3Array> Vgmc = dynamic_cast<osg::Vec3Array *> (gmC->getVertexArray());
	osg::Vec3 firstVertex = Vgmc->front();
	int numVers = Vgmc->size();
	int iFlag = 0;
	for(std::vector<LakeEdgeItem>::iterator p = allEdgedLakes.begin(); p != allEdgedLakes.end(); ++p)
	{
		if((numVers == p->numOfVertex) && (firstVertex == p->firstVertex))
		{
			iFlag = 1; break;
		}
	}
	if(!iFlag)
	{
		getEdgesOfGeometry(&allLakeEdges, gmC);
		LakeEdgeItem lei;	lei.numOfVertex = numVers; lei.firstVertex = firstVertex; allEdgedLakes.push_back(lei);
	}
	/* Add by Galen 2013-02-22 End */
	
	
	
	/* 检查地形面的每个顶点。 */
	for( std::vector<FindValidVertex>::iterator v = faceVertice.begin(); v != faceVertice.end(); ++v )
	{
		double minDistance = 999999.999;
		std::vector<osg::Vec3> minDistVxs;
		
		/* 对于地形面上的每个点，首先要不与gmC相交 */
		osg::Vec3 isect;
		if( TerrainProcessingUtils::getIntersect( flatFeatureGd, v->pVF, isect ))	
		{
			/* Add by Galen 2013-02-25 把这个点的高程直接赋予湖面高程，不管后面是否有用过。*/
			Geodesy geo;
			Carton  cat;
			osg::Vec3 vx = v->pVF;

			cat.x = vx.x(); cat.y = vx.y(); cat.z = vx.z();
			CartToGeod(&cat, &geo);	geo._elevation = elevationOfLake; GeodToCart(&geo, &cat);
			vx.x() = cat.x; vx.y() = cat.y; vx.z() = cat.z;
			v->pV = vx;
			terrainVertices[v->idx].pV = v->pV;
			/* Add by Galen 2013-02-25 End */
			continue;
		}

		/* 计算边集合kept 里面所有的边的两个端点到 *v 点的距离之和的最小值，这个值最小说明这条边距离该点最近。*/
		for( std::vector<TerrainProcessingUtils::Edge>::iterator p = allLakeEdges.begin(); p != allLakeEdges.end(); p++ )
		{
			double dist = getDistanceFromTwoVertexToVertex(p->A, p->B, v->pVF);
			if(dist < minDistance)
			{
				minDistance = dist;	minDistVxs.clear();	minDistVxs.push_back(p->A);	minDistVxs.push_back(p->B);
			}
		}
		
		/* 下面计算边 minDistVxs， 与点 *v 之间的距离。*/
		/* 面积高度检验法， */
		/* 首先求出这个三角形的面积 */
		double area = getArea(minDistVxs[0], minDistVxs[1], v->pVF);
		/* 求出文化信息边的长度 */
		double maxLen = distanceV(minDistVxs[0], minDistVxs[1]);
		double heigh = 2 * area / maxLen;
		
		/* 记录这个点的信息 */
		if(heigh < MINDISTANCETOEDGE)
		{
			VertexDistance vd;	vd.pVF = v->pVF;	vd.distance = heigh;	vertexNearLake.push_back(vd);
		}
	}	// v

}
/* Add by Galen 2013-02-20  End */


/* 重构features, 把 features 上不跟地景相交的面删掉。
   目前这个算法运算量太大了。
   可以重新写一个算法：
   首先计算地景面的包围球，
   再计算每个三角面的包围球，
   只要两个球心之间的距离大于两球半径之和，就说明面不相交。
   如果两个球心之间的距离小于两球半径之和，那么再用插值的办法来仔细判断。
 */
bool reConstructFeature(osg::Geode* gd, osg::Geometry *g, osg::ref_ptr<osg::Vec3Array> tbLoop)
{
	/* 首先判断几何体g整体的包围球是否与地形gd包围球有重合 */
	const osg::BoundingSphere bsTerrain  = gd->getBound();
	const osg::BoundingSphere bsLineloop = g->getBound();

	/* 计算两个包围球球心之间的距离，和两个包围球半径之和 */	
	float sumRadius = bsTerrain.radius() + bsLineloop.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsLineloop.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return false;

	/* 再判断几何体g内每个三角形的包围球是否与地形gd包围球有重合*/
	osg::ref_ptr<osg::Vec3Array> coords = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	osg::ref_ptr<osg::DrawElementsUInt> deui = new osg::DrawElementsUInt(GL_TRIANGLES);
	bool bvalid = true;
	for (unsigned int ipr=0; ipr< g->getNumPrimitiveSets(); ipr++) 
	{
		osg::ref_ptr<osg::PrimitiveSet> prset = g->getPrimitiveSet(ipr);

		if(prset->getMode() == osg::PrimitiveSet::LINE_LOOP)
		{
		}
		else
		{
			switch( prset->getType() )
			{
				case osg::PrimitiveSet::DrawArraysPrimitiveType:
					{
						osg::DrawArrays* gdui = dynamic_cast<osg::DrawArrays*>(prset.get());
						if( gdui )
						{
							switch(prset->getMode() )
							{
								case osg::PrimitiveSet::TRIANGLES:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
									{
										unsigned int a,b,c;
										a = gdui->index(i);
										b = gdui->index(i+1);
										c = gdui->index(i+2);	

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
									}
									break;
								case osg::PrimitiveSet::QUADS:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=4)
									{
										unsigned int a,b,c,d;
										a = gdui->index(i);
										b = gdui->index(i+1);
										c = gdui->index(i+2);
										d = gdui->index(i+3);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}

										if (checkTriwithTerrain(*(coords->begin()+c), *(coords->begin()+d), *(coords->begin()+a), gd, tbLoop)!=0)
										{
											deui->push_back(c);
											deui->push_back(d);
											deui->push_back(a);
										}
									}
									break;
								default:
									break;
							}
						}
					 }
					break;
				case osg::PrimitiveSet::DrawArrayLengthsPrimitiveType:
					break;
				case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
					{
						osg::DrawElementsUByte *gdui = dynamic_cast<osg::DrawElementsUByte *>(prset.get());
						unsigned int a,b,c;
						if( gdui )
						{
							switch(prset->getMode() )
							{
								case osg::PrimitiveSet::TRIANGLES:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
									{
										a = gdui->getElement(i);
										b = gdui->getElement(i+1);
										c = gdui->getElement(i+2);		

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_STRIP:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										a=b; b=c;
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_FAN:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										b=c;
									}
									break;
								case osg::PrimitiveSet::QUADS:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=4)
									{
										unsigned int a,b,c,d;
										a = gdui->index(i);
										b = gdui->index(i+1);
										c = gdui->index(i+2);
										d = gdui->index(i+3);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}

										if (checkTriwithTerrain(*(coords->begin()+c), *(coords->begin()+d), *(coords->begin()+a), gd, tbLoop)!=0)
										{
											deui->push_back(c);
											deui->push_back(d);
											deui->push_back(a);
										}
									}
									break;
								default:
									break;
							}
						}
					}
					break;
				case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
					{
						osg::DrawElementsUShort *gdui = dynamic_cast<osg::DrawElementsUShort *>(prset.get());
						unsigned int a,b,c;
						if( gdui )
						{
							switch(prset->getMode() )
							{
								case osg::PrimitiveSet::TRIANGLES:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
									{
										a = gdui->getElement(i);
										b = gdui->getElement(i+1);
										c = gdui->getElement(i+2);		

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_STRIP:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										a=b; b=c;
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_FAN:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										b=c;
									}
									break;
								case osg::PrimitiveSet::QUADS:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=4)
									{
										unsigned int a,b,c,d;
										a = gdui->index(i);
										b = gdui->index(i+1);
										c = gdui->index(i+2);
										d = gdui->index(i+3);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}

										if (checkTriwithTerrain(*(coords->begin()+c), *(coords->begin()+d), *(coords->begin()+a), gd, tbLoop)!=0)
										{
											deui->push_back(c);
											deui->push_back(d);
											deui->push_back(a);
										}
									}
									break;
								default:
									break;
							}
						}
					}
					break;
				case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
					{
						osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
						unsigned int a,b,c;
						if( gdui )
						{
							switch(prset->getMode() )
							{
								case osg::PrimitiveSet::TRIANGLES:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
									{
										a = gdui->getElement(i);
										b = gdui->getElement(i+1);
										c = gdui->getElement(i+2);		

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_STRIP:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										a=b; b=c;
									}
									break;
								case osg::PrimitiveSet::TRIANGLE_FAN:
									a = gdui->getElement(0);
									b = gdui->getElement(1);
									for(unsigned int i=2; i<gdui->getNumIndices(); i++)
									{
										c = gdui->getElement(i);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}
										b=c;
									}
									break;
								case osg::PrimitiveSet::QUADS:
									for(unsigned int i=0; i<gdui->getNumIndices(); i+=4)
									{
										unsigned int a,b,c,d;
										a = gdui->index(i);
										b = gdui->index(i+1);
										c = gdui->index(i+2);
										d = gdui->index(i+3);

										if (checkTriwithTerrain(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), gd, tbLoop)!=0)
										{
											deui->push_back(a);
											deui->push_back(b);
											deui->push_back(c);
										}

										if (checkTriwithTerrain(*(coords->begin()+c), *(coords->begin()+d), *(coords->begin()+a), gd, tbLoop)!=0)
										{
											deui->push_back(c);
											deui->push_back(d);
											deui->push_back(a);
										}
									}
									break;
								default:
									break;
							}
						}
					}
					break;
				default:
					break;
			}	//case
		}
	}	//for

	if(!(deui->empty()))
	{
		g->removePrimitiveSet(0, g->getNumPrimitiveSets());
		g->addPrimitiveSet(deui);
	}else
		bvalid = false;

	return bvalid;
}

/* 计算一个多边形(线环)的面积。*/
double GetAreaOfLineloop(osg::ref_ptr<osg::Vec3Array> ll)
{
	/* 先把当前perlist线环都三角化成一个面，为了跟边界段集合做相交检查 */				
	osg::ref_ptr<osg::Geometry> g = ConvertFromLINELOOPToTranglesFace(ll);
	if(g == NULL)	return 0.0;
	
	/* 计算所有的三角面的面积之和。*/
	double totalArea = 0.0;
	osg::ref_ptr<osg::Vec3Array> coords = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	for (unsigned int ipr=0; ipr< g->getNumPrimitiveSets(); ipr++) 
	{
		osg::ref_ptr<osg::PrimitiveSet> prset = g->getPrimitiveSet(ipr);
		switch( prset->getType() )
		{
			case osg::PrimitiveSet::DrawArraysPrimitiveType:
			case osg::PrimitiveSet::DrawArrayLengthsPrimitiveType:
			case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
			case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
				break;
			case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
				{
					osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
					if( gdui )
					{
						switch(prset->getMode() )
						{
							case osg::PrimitiveSet::TRIANGLES:
								for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
								{
									unsigned int a,b,c;
									a = gdui->getElement(i);
									b = gdui->getElement(i+1);
									c = gdui->getElement(i+2);
									double area = getArea(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c));
									totalArea += area;
								}
								break;
							default:
								break;
						}
					}
				}
				break;
		}	//case
	}	//for
	return totalArea;
}


/* 把一个包含有多个primetives的dc合并到dclist里面 */
void MergeDcIntoDclist(std::vector<osg::ref_ptr<ArealConstraint> > *dclist, osg::ref_ptr<ArealConstraint> mergedc, string dcname)
{
	/* 先把它合并成一个gmlist */
	osg::ref_ptr<osgUtil::DelaunayConstraint> dc = mergedc;
	osg::ref_ptr<osg::Vec3Array> coords = dynamic_cast<osg::Vec3Array *>(dc->getVertexArray());

	std::vector<osg::ref_ptr<osg::Vec3Array>> perList;
	for(unsigned int i=0; i<mergedc->getNumPrimitiveSets(); i++)
	{
		osg::ref_ptr<osg::PrimitiveSet> prset = mergedc->getPrimitiveSet(i);

		if(prset->getMode() == osg::PrimitiveSet::LINE_LOOP)
		{
			osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
			switch( prset->getType() )
			{
			case osg::PrimitiveSet::DrawArraysPrimitiveType:
				break;
			case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
				{
					osg::DrawElementsUByte *gdui = dynamic_cast<osg::DrawElementsUByte *>(prset.get());
					if( gdui )
					{
						for(unsigned int i=0; i<gdui->getNumIndices(); i++)
						{
							unsigned char a = gdui->getElement(i);
							vx->push_back(*(coords->begin()+a));
						}
						
					}
				}
				break;
			case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
				{
					osg::DrawElementsUShort *gdui = dynamic_cast<osg::DrawElementsUShort *>(prset.get());
					if( gdui )
					{
						for(unsigned int i=0; i<gdui->getNumIndices(); i++)
						{
							unsigned short a = gdui->getElement(i);
							vx->push_back(*(coords->begin()+a));
						}
						
					}
				}
				break;
			case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
				{
					osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
					if( gdui )
					{
						for(unsigned int i=0; i<gdui->getNumIndices(); i++)
						{
							unsigned int a = gdui->getElement(i);
							vx->push_back(*(coords->begin()+a));
						}
						
					}
				}
				break;
			}

			perList.push_back(vx);
		}

	}

	/* 如果是湖泊，那么计算每片湖泊的面积，如果某个湖面面积很小很小，这个湖面可忽略不计 */
	if(dcname == "Lake")
	{
		std::vector<osg::ref_ptr<osg::Vec3Array>> temp;
		for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perList.begin(); q != perList.end(); q++ )
			temp.push_back( *q );
		perList.clear();
		for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= temp.begin(); q != temp.end(); q++ )
		{
			double area = GetAreaOfLineloop(*q);
			if(area > AREA_LIMIT)	perList.push_back(*q);
		}
	}
	
	/* 把线环变成dc压入到dclist */
	for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perList.begin(); q != perList.end(); q++ )
	{
		if (q->get()->size()>=3) 
		{
			osg::ref_ptr<ArealConstraint> dc = new ArealConstraint;
			dc.get()->setName(dcname);
			dc->setVertexArray( q->get() );
			dc->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,q->get()->size()) );
			dclist->push_back( dc.get() );
		}
	}	//PerList

}


/* 把所有的具有相同名称的dc合并成为一个，并计算覆盖情况 */
void MergeDclistStepOne(std::vector<osg::ref_ptr<ArealConstraint> > *dclist)
{
	if(dclist->size() == 0)	return;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* 下面这段代码针对 文化信息路有交叉的情况的处理。*/
	/* 先统计有多少个同名的 dc，多于1个以上的，才做 merge */
	totalRoad = 0;
	totalIsland = 0;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
	{
		if (p->get()->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING) totalRoad++;
		if (p->get()->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) totalIsland++;
	}

	/* 把所有的 名称为 Road 的合并成一个 Road */
	if(totalRoad > 1)
	{
		osg::ref_ptr<ArealConstraint> mergedc;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
		{
			if (p->get()->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)
			{
				if (mergedc.valid()==false)
					mergedc = new ArealConstraint(*(p->get()));
				else
					mergedc->merge(p->get());
			}
		}
		if (mergedc.valid())
		{
			mergedc->handleOverlaps();
			mergedc.get()->setName("Road");
		}
	
		/* 从mergedc里面找出新增的顶点。*/
		osg::ref_ptr<osg::Vec3Array> dcVry = dynamic_cast<osg::Vec3Array *> (mergedc->getVertexArray());
		for( osg::Vec3Array::iterator r = dcVry->begin(); r != dcVry->end(); r++ )
			addVertexIntoTotalFirstSearth(*r);


		osg::ref_ptr<ArealConstraint> mergedc2;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
		{
			if (p->get()->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)
			{
				if (mergedc2.valid()==false)
					mergedc2 = new ArealConstraint(*(p->get()));
				else
					mergedc2->merge(p->get());
			}
		}

		if (mergedc2.valid())
		{
			mergedc2->handleOverlaps1();
			mergedc2.get()->setName("Road");
		}
	
		if (mergedc.valid() || mergedc2.valid())
		{
			std::vector<osg::ref_ptr<ArealConstraint> > tmp_dclist;
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
			{
				if (p->get()->getName().substr(0, strlen(DEFINE_ROAD_STRING)) != DEFINE_ROAD_STRING)
					tmp_dclist.push_back(p->get());
			}

			dclist->clear();
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= tmp_dclist.begin(); p != tmp_dclist.end(); p++ )
				dclist->push_back(p->get());
			
			/* 合并mergedc */
			if(mergedc.valid())		MergeDcIntoDclist(dclist, mergedc, "Road");
			if(mergedc2.valid())	MergeDcIntoDclist(dclist, mergedc2, "Road");
		}
	}	//totalRoad > 1;
}


/* 把Lake和Island的dc合并成为一个，并计算覆盖情况 */
void MergeLakeAndIsland(std::vector<osg::ref_ptr<ArealConstraint> > *dclist)
{
	if(dcLakes.size() == 0) return;
	
	/* 先统计有多少个"Island" 或者 “Lake” dc，多于1个以上的，才做 merge */
	totalLake = dcLakes.size();
	totalIsland = 0;
	if(dclist->size() > 0)
	{
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
		{
			if (p->get()->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) totalIsland++;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* 下面这段代码针对 文化信息 大海中有若干岛屿 的情况的处理。*/
	/* 把要挖的文化信息river和 不要挖的 island 之间也做 merge*/
	if(((totalLake > 0) && (totalIsland > 0)) || (totalLake > 1))
	{
		osg::ref_ptr<ArealConstraint> mergedc;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
		{
			if (mergedc.valid()==false)
				mergedc = new ArealConstraint(*(p->get()));
			else
				mergedc->merge(p->get());
		}
		
		if(totalIsland > 0)
		{
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
			{
				if (p->get()->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING)
				{
					mergedc->merge(p->get());
				}
			}
		}
		
		if (mergedc.valid())
		{
			mergedc->handleOverlaps();
			mergedc.get()->setName("Lake");
			
			if(mergedc.valid())
			{
				if(mergedc->getNumPrimitiveSets() != 0)
				{
					dcLakes.clear();
					MergeDcIntoDclist(&dcLakes, mergedc, "Lake");
				}
			}
		}
	}
}



/* 创建同动态草的面相交的点的序列 */
void CreateValidVertexOnFaceOfDynamicGrass(MaxMinLonLat *mmll, osg::ref_ptr<osg::Geometry> gm, double elevMid, double elevOfGrass)
{
	/* 计算出 1 米所占据的经纬度是多少 deltaLon, deltaLat */
	Geodesy p1, p2, geo;
	Carton car;
	p1._lon = mmll->minLong; p1._lat = mmll->maxLati; p1._elevation = 0;
	p2._lon = mmll->maxLong; p2._lat = mmll->maxLati; p2._elevation = 0;
	double distLon = distanceM(&p1, &p2);
	p1._lon = mmll->minLong; p1._lat = mmll->minLati; p1._elevation = 0;
	p2._lon = mmll->minLong; p2._lat = mmll->maxLati; p2._elevation = 0;
	double distLat = distanceM(&p1, &p2);
	double deltaLon = (mmll->maxLong - mmll->minLong) / distLon;
	double deltaLat = (mmll->maxLati - mmll->minLati) / distLat;
	osg::ref_ptr<osg::Geode> gFace = new osg::Geode;
	gFace->addDrawable(gm);

	/* 创建动态草结点 */
	osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
	gnode->setName("DynamicGrass");
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
	gnode->addDrawable(geom);
	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	geom->setVertexArray(vx);

	/* 插入点 */
	for(double i = mmll->minLong; i < mmll->maxLong; i += UNITDISTANCE * deltaLon)
	{
		for(double j = mmll->minLati; j < mmll->maxLati; j += UNITDISTANCE * deltaLat)
		{
			/* 经纬高转换成xyz */
			geo._lon = i; geo._lat = j; geo._elevation = elevMid;	GeodToCart(&geo, &car);
			osg::Vec3 currVert(car.x, car.y, car.z);

			/* 计算每个点是否跟 g 面相交，求交点 */
			osg::Vec3 isect;
			if( TerrainProcessingUtils::getIntersect( gFace, currVert, isect ))
			{
				if(elevOfGrass == 0.0)
					vx->push_back(isect);
				else
				{
					car.x = isect.x();	car.y = isect.y(); car.z = isect.z();
					CartToGeod(&car, &geo);	geo._elevation = elevOfGrass;
					GeodToCart(&geo, &car);
					osg::Vec3 cVert(car.x, car.y, car.z);
					vx->push_back(cVert);
				}
			}
		}
	}
	
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, vx->size()) );
	dynamicGrass_gdlist.push_back(gnode); 

	return;
}

/* 创建动态草的面 */
void CreateFaceOfDynamicGrass(osgTDS::TexCoordPreserver *tcp)
{
	if(dclist.empty()) return;

	//创建一个新的点集合，只包含面上的点
	osg::ref_ptr<osg::Vec3Array> geVx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		geVx->push_back(f->pVF);

	osg::ref_ptr<osg::Vec3Array> dcoords = TerrainProcessingUtils::uniqueCoords(geVx.get());
	osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(dcoords.get());
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		if(p->get()->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) == DEFINE_DYNAMICGRASS_STRING)
			dt->addInputConstraint( p->get() );
	}

	dt->triangulate();
	
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		if(p->get()->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) == DEFINE_DYNAMICGRASS_STRING)
			dt->removeInternalTriangles( p->get() );
	}


	double elevOfGrass = 0.0;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		std::string name = p->get()->getName();
		if(p->get()->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) != DEFINE_DYNAMICGRASS_STRING)	continue;
		
		/* 找出当前动态草的高度 */
		SplitName(name);
		if(NameParts.size() >= 2)	elevOfGrass = atof(NameParts[1].c_str());

		if (p->get()->getinteriorTrisSize()>0)
		{
			if(1)
			{
				osg::ref_ptr<osg::Vec3Array> arpts=p->get()->getPoints(dcoords.get());

				/* 把三角化后的所有点集合都存在validVertices里面 */
				/* 并把validVertices里面所有点，由投影点映射到空间点 */
				std::vector<ProjectionMap> validProjection;
				MaxMinLonLat mmll;
				InitMMLLEx(&mmll, "");
				double elevTotal = 0.0;

				int j = 0;
				validVertices.clear();	validElement.clear();
				for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
				{
					FindValidVertex vVert;
					osg::Vec3 pSpace = getSpaceVertex(*r);
					vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
					validVertices.push_back(vVert);
					
					ProjectionMap pm;
					getProjectionMapByVertex(&pm, pSpace);
					SetMMLLEx(&mmll, pm.lon, pm.lat);	elevTotal += pm.elev;
				}
				double elevMid = elevTotal / j;

				//找出所有的有效面。存放在 validElement 里面
				osg::DrawElementsUInt *validdui = p->get()->getTriangles();
				for(unsigned int i=0; i<validdui->getNumIndices(); i+=3)
				{
					unsigned int a,b,c;
					a = validdui->getElement(i);
					b = validdui->getElement(i+1);
					c = validdui->getElement(i+2);
					if(	(checkTriOnBorder(&validVertices[a], &validVertices[b], &validVertices[c])) &&
							(checkTriOnDclist(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF)) )
					{
						osg::Vec3 normal = getNormal(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF);
						if(normal.z() > 0.0)
						{
							validElement.push_back(a);
							validElement.push_back(b);
							validElement.push_back(c);
						}else{
							validElement.push_back(a);
							validElement.push_back(c);
							validElement.push_back(b);
						}
					}
				}

				//创建一个新的点阵列，加入所有有效面
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//把只包含合法三角形的validElement里面的元素加进来
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

				/* 创建最终完成切割后的几何面。是全空间点立体面。*/
				osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
				gm->setVertexArray( validCoords.get() );
				gm->addPrimitiveSet( dui );;

				/* 找出所有与切割后几何面相交的点。点与点的间距是2米*/
				/* 根据计算出来的动态草的面，计算出内部的插入点 */
				CreateValidVertexOnFaceOfDynamicGrass(&mmll, gm, elevMid, elevOfGrass);
				
			}
		}
	} 
}



/* 把细节纹理的dc单独拿出来。 */
void MoveOutDetailTxtDc()
{
	/* 先把湖面拿出来。*/
	std::vector<osg::ref_ptr<ArealConstraint> > noLakes;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		std::string name= p->get()->getName();
		if (name.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			dcLakes.push_back(*p);
		else
			noLakes.push_back(*p);
	}
	
	dclist.clear();
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= noLakes.begin(); p != noLakes.end(); p++ )
		dclist.push_back(*p);
	
	
	/* 再把细节纹理面拿出来 */	
	std::vector<osg::ref_ptr<ArealConstraint> > nodtl;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		std::string name= p->get()->getName();
		if (isDetailedtxt(name)==true)
			dcDetail.push_back(*p);
		else
			nodtl.push_back(*p);
	}
	
	dclist.clear();
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= nodtl.begin(); p != nodtl.end(); p++ )
		dclist.push_back(*p);
}


/* 创建湖的面 */
void CreateFaceOfLake(osgTDS::TexCoordPreserver *tcp)
{
	if(dcLakes.empty()) return;


	//创建一个新的点集合，只包含面上的点
	osg::ref_ptr<osg::Vec3Array> geVx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		geVx->push_back(f->pVF);

	osg::ref_ptr<osg::Vec3Array> dcoords = TerrainProcessingUtils::uniqueCoords(geVx.get());
	osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(dcoords.get());
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
		dt->addInputConstraint( p->get() );

	dt->triangulate();
	
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
	{
		dt->removeInternalTriangles( p->get() );
	} 


	osg::ref_ptr<osg::Geode> gnode = NULL;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
	{
		std::string name= p->get()->getName();
		
		/* 如果name行尾有空格，则把它清除掉 */
		char namestr[128];
		sprintf(namestr, "%s", name.c_str());
		for(unsigned int i=0; i<strlen(namestr); i++)
		{
			if((namestr[i]==' ') || (namestr[i]==0x9))  namestr[i] = 0;
		}
		name = namestr;

		if (p->get()->getinteriorTrisSize()>0)
		{
			if (gnode==NULL || gnode->getName()!=name)
			{
				gnode = new osg::Geode;
				osg::ref_ptr<osg::StateSet> st = new osg::StateSet(*(pCurTerrain->getStateSet()));
				osg::Texture2D* t2d = dynamic_cast<osg::Texture2D*>(st->getTextureAttribute(0,osg::StateAttribute::TEXTURE));
				st->setTextureAttributeAndModes(0,t2d,osg::StateAttribute::ON);
				gnode->setStateSet(st.get());
				gnode->setName(name);
				detailed_gdlist.push_back(gnode); 
			}

			if(1)
			{
				/* 找出有效顶点集合 */
				osg::ref_ptr<osg::Vec3Array> arpts=p->get()->getPoints(dcoords.get());

				/* 把三角化后的所有点集合都存在validVertices里面 */
				/* 并把validVertices里面所有点，由投影点映射到空间点 */
				int j = 0;
				validVertices.clear();	validElement.clear();
				for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
				{
					FindValidVertex vVert;
					osg::Vec3 pSpace = getSpaceVertex(*r);
					vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
					validVertices.push_back(vVert);
				}

				//找出所有的有效面。存放在 validElement 里面
				osg::DrawElementsUInt *validdui = p->get()->getTriangles();
				for(unsigned int i=0; i<validdui->getNumIndices(); i+=3)
				{
					unsigned int a,b,c;
					a = validdui->getElement(i);
					b = validdui->getElement(i+1);
					c = validdui->getElement(i+2);
					if(	(checkTriOnBorder(&validVertices[a], &validVertices[b], &validVertices[c])) &&
							(checkTriOnDclist(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF)) )
					{
						osg::Vec3 normal = getNormal(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF);
						if(normal.z() > 0.0)
						{
							validElement.push_back(a);
							validElement.push_back(b);
							validElement.push_back(c);
						}else{
							validElement.push_back(a);
							validElement.push_back(c);
							validElement.push_back(b);
						}
					}
				}

				//创建一个新的点阵列，加入所有有效面
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//把只包含合法三角形的validElement里面的元素加进来
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

#ifdef ADD_SKIRT
				/* 创建裙边 */
				CreateSkirt(validCoords, dui, 0);
#endif

				/* 创建最终完成切割后的几何体。*/
				osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
				gm->setVertexArray( validCoords.get() );
				osg::ref_ptr<osg::Vec2Array> getc = tcp->getTexCoords( *(validCoords.get()) );
				gm->setTexCoordArray( 0, getc.get() );
				// Restore normal 
				osg::ref_ptr<osg::Vec3Array> nc = tcp->getNormals(  *(validCoords.get()) );
				gm->setNormalArray( nc.get() );
				gm->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
				// Add an overall white color
				osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
				color->push_back(osg::Vec4(1,1,1,1));
				gm->setColorArray(color.get());
				gm->setColorBinding( osg::Geometry::BIND_OVERALL );
				gm->addPrimitiveSet( dui );;

#ifdef ADD_SKIRT
				for( std::vector< osg::ref_ptr<osg::DrawElementsUInt> >::iterator p = cultSkirtDrawElements.begin(); p != cultSkirtDrawElements.end(); p++ )
					gm->addPrimitiveSet( *p );
#endif

				gnode->addDrawable(gm.get());
			}
		}
	} 
	
	
}

/* 创建细节纹理的面，无文化信息的 */
void CreateFaceOfDetailTxtWhenNoCulture(osgTDS::TexCoordPreserver *tcp)
{
	if(dcDetail.empty()) return;

	//创建一个新的点集合，只包含面上的点
	osg::ref_ptr<osg::Vec3Array> geVx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		geVx->push_back(f->pVF);

	osg::ref_ptr<osg::Vec3Array> dcoords = TerrainProcessingUtils::uniqueCoords(geVx.get());
	osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(dcoords.get());
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
		dt->addInputConstraint( p->get() );

	dt->triangulate();
	
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
	{
		dt->removeInternalTriangles( p->get() );
	} 


	osg::ref_ptr<osg::Geode> gnode = NULL;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
	{
		std::string name= p->get()->getName();
		
		/* 如果name行尾有空格，则把它清除掉 */
		char namestr[128];
		sprintf(namestr, "%s", name.c_str());
		for(unsigned int i=0; i<strlen(namestr); i++)
		{
			if((namestr[i]==' ') || (namestr[i]==0x9))  namestr[i] = 0;
		}
		name = namestr;

		if (p->get()->getinteriorTrisSize()>0)
		{
			if (gnode==NULL || gnode->getName()!=name)
			{
				gnode = new osg::Geode;
				osg::ref_ptr<osg::StateSet> st = new osg::StateSet(*(pCurTerrain->getStateSet()));
				osg::Texture2D* t2d = dynamic_cast<osg::Texture2D*>(st->getTextureAttribute(0,osg::StateAttribute::TEXTURE));
				st->setTextureAttributeAndModes(0,t2d,osg::StateAttribute::ON);
				gnode->setStateSet(st.get());
				gnode->setName(name);
				detailed_gdlist.push_back(gnode); 
			}

			if(1)
			{
				osg::ref_ptr<osg::Vec3Array> arpts=p->get()->getPoints(dcoords.get());

				/* 把三角化后的所有点集合都存在validVertices里面 */
				/* 并把validVertices里面所有点，由投影点映射到空间点 */
				int j = 0;
				validVertices.clear();	validElement.clear();
				for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
				{
					FindValidVertex vVert;
					osg::Vec3 pSpace = getSpaceVertex(*r);
					vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
					validVertices.push_back(vVert);
				}

				//找出所有的有效面。存放在 validElement 里面
				osg::DrawElementsUInt *validdui = p->get()->getTriangles();
				for(unsigned int i=0; i<validdui->getNumIndices(); i+=3)
				{
					unsigned int a,b,c;
					a = validdui->getElement(i);
					b = validdui->getElement(i+1);
					c = validdui->getElement(i+2);
					if(	(checkTriOnBorder(&validVertices[a], &validVertices[b], &validVertices[c])) &&
							(checkTriOnDclist(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF)) )
					{
						osg::Vec3 normal = getNormal(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF);
						if(normal.z() > 0.0)
						{
							validElement.push_back(a);
							validElement.push_back(b);
							validElement.push_back(c);
						}else{
							validElement.push_back(a);
							validElement.push_back(c);
							validElement.push_back(b);
						}
					}
				}

				//创建一个新的点阵列，加入所有有效面
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//把只包含合法三角形的validElement里面的元素加进来
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

#ifdef ADD_SKIRT
				/* 创建裙边 */
				CreateSkirt(validCoords, dui, 0);
#endif

				/* 创建最终完成切割后的几何体。*/
				osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
				gm->setVertexArray( validCoords.get() );
				osg::ref_ptr<osg::Vec2Array> getc = tcp->getTexCoords( *(validCoords.get()) );
				gm->setTexCoordArray( 0, getc.get() );
				// Restore normal 
				osg::ref_ptr<osg::Vec3Array> nc = tcp->getNormals(  *(validCoords.get()) );
				gm->setNormalArray( nc.get() );
				gm->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
				// Add an overall white color
				osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
				color->push_back(osg::Vec4(1,1,1,1));
				gm->setColorArray(color.get());
				gm->setColorBinding( osg::Geometry::BIND_OVERALL );
				gm->addPrimitiveSet( dui );;

#ifdef ADD_SKIRT
				for( std::vector< osg::ref_ptr<osg::DrawElementsUInt> >::iterator p = cultSkirtDrawElements.begin(); p != cultSkirtDrawElements.end(); p++ )
					gm->addPrimitiveSet( *p );
#endif

				gnode->addDrawable(gm.get());
			}
		}
	} 
}

/* 创建细节纹理的面，有文化信息的 */
void CreateFaceOfDetailTxtWhenHasCulture(osg::ref_ptr<osg::Geometry> geom)
{
	if(dcDetail.empty()) return;

	osg::ref_ptr<osg::Geode> gnode = NULL;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
	{
		std::string name= p->get()->getName();
		
		/* 如果name行尾有空格，则把它清除掉 */
		char namestr[128];
		sprintf(namestr, "%s", name.c_str());
		for(unsigned int i=0; i<strlen(namestr); i++)
		{
			if((namestr[i]==' ') || (namestr[i]==0x9))  namestr[i] = 0;
		}
		name = namestr;

		if (gnode==NULL || gnode->getName()!=name)
		{
			gnode = new osg::Geode;
			osg::ref_ptr<osg::StateSet> st = new osg::StateSet(*(pCurTerrain->getStateSet()));
			osg::Texture2D* t2d = dynamic_cast<osg::Texture2D*>(st->getTextureAttribute(0,osg::StateAttribute::TEXTURE));
			st->setTextureAttributeAndModes(0,t2d,osg::StateAttribute::ON);
			gnode->setStateSet(st.get());
			gnode->setName(name);
			detailed_gdlist.push_back(gnode); 
		}

		if(1)
		{
			osg::ref_ptr<osg::Geometry> gm = new osg::Geometry(*geom, osg::CopyOp::DEEP_COPY_ALL);
			gnode->addDrawable(gm.get());
		}
	} 
}

/* 检查dclist情况，是不是湖泊和岛屿具有相同的dc。满足以下要求： 
   1、dclist只有两个dc，
   2、两个dc顶点序列是一样的。
   3、两个dc的名字也一样。
   动作：舍弃其中一个dc，只留一个。
*/
int CheckDclistAndMerge(std::vector<osg::ref_ptr<ArealConstraint> > *dclist)
{
	int ret = 0;
	if(dclist->size() != 2) return ret;
	osg::ref_ptr<osgUtil::DelaunayConstraint> dc0 = dclist->at(0);
	osg::ref_ptr<osgUtil::DelaunayConstraint> dc1 = dclist->at(1);
	osg::Vec3Array *dc0coords = dynamic_cast<osg::Vec3Array *>(dc0->getVertexArray());
	osg::Vec3Array *dc1coords = dynamic_cast<osg::Vec3Array *>(dc1->getVertexArray());
	if(dc0coords->getNumElements() != dc1coords->getNumElements() ) return ret;

	osg::Vec3Array::iterator iter_dc0 = dc0coords->begin();
	osg::Vec3Array::iterator iter_dc1 = dc1coords->begin();
	for(; iter_dc0 != dc0coords->end() && iter_dc1 != dc1coords->end() ; iter_dc0++)
	{
		osg::Vec3 v0 = *iter_dc0;
		osg::Vec3 v1 = *iter_dc1;	iter_dc1++;
		if(v0 != v1) return ret;
	}
	
	/* 如果顶点序列相同，名字也相同 */
	if(dc0->getName() == dc1->getName())
	{
		dclist->clear();
		osg::ref_ptr<ArealConstraint> dc = dynamic_cast<ArealConstraint *> (dc0.get());
		dclist->push_back(dc);
	}

	return ret;
}

/* 初始化所有用到的容器 */
void InitVectors()
{
	totalVertices.clear();
	terrainVertices.clear();
	terrainBorderVertices.clear();
	terrainInnerVertices.clear();
	allBorderVertices.clear();
	faceVertice.clear();
	skirtVertics.clear();
	faceElement.clear();
	skirtElement.clear();
	intersectlist.clear();
	roadTexCoordSave.clear();
	validVertices.clear();
	validBorderVertices.clear();
	checkBorderVertices.clear();
	invalidVertices.clear();
	innerVerticesInterWithCulture.clear();
	validElement.clear();
	perimeterList.clear();
	validCSM.clear();
	dclist.clear();
	dcDetail.clear();
	dcLakes.clear();
	if(IsTopTerrain)	IntersectVertices.clear();
	vertexNearLake.clear();
}

/* 将整个地景面上的所有点分为两个部分，一部分是地景面上的点，一部分是裙边上的点。*/
/* 根据纹理坐标情况分发顶点信息到faceVertices和skirtVertices。分发原则是，与上一个点纹理信息相同但高度低于上一个点的，就是skirt点。 */
void SplitTerrainVertexToFaceAndSkirt()
{
	/* 首先第一个点肯定是face面上的点。 */
	faceVertice.push_back(terrainVertices.front());
	
	/* 再按照纹理坐标相同，顶点高度高低来区分是 裙边点 还是 面上点 */
	for( std::vector<FindValidVertex>::iterator p= terrainVertices.begin()+1; p != terrainVertices.end(); p++ )
	{
		int isSkirt = 0;
		for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		{
			if((p->pT == f->pT) && (p->pV[2] < f->pV[2]))
			{
				skirtVertics.push_back(*p);
				isSkirt = 1;
			}
			if((p->pT == f->pT) && (p->pV[2] > f->pV[2]))
			{
				skirtVertics.push_back(*f);
				*f = *p;
				isSkirt = 1;
			}
		}
		if(!isSkirt)
			faceVertice.push_back(*p);
	}
}

/* 将边界线环上的所有点，由裙边点，都替换成面上的点 */
void ChangeBorderVertexFromSkirtToFace(osg::Geometry *g)
{
	std::vector<osg::ref_ptr<osg::Vec3Array> > periList;

	/* 计算出地景边缘的线环。这些点有可能都是裙边上的点。*/
	TerrainProcessingUtils::computePerimeterAll(g, periList);
	if(periList.size() > 0) terrainBorderLineLoop = periList[0];


	/* 将border线环上的所有顶点信息读入到 terrainBorderVertices 数组。*/
	unsigned int j = 0;
	for( osg::Vec3Array::iterator r = terrainBorderLineLoop->begin(); r != terrainBorderLineLoop->end(); r++ )
	{
		FindValidVertex fvv;
		fvv.pV = *r; fvv.idx = j; fvv.newIdx = 0; j++; fvv.isValid = -1;
		terrainBorderVertices.push_back(fvv);
	}
    
	/* 根据 terrainBorderVertices 各个顶点的pV值 和 skirtVertics 的pV值一一对应，找出 terrainBorderVertices 里面对应的 pT值 。*/
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
	{
		for( std::vector<FindValidVertex>::iterator f = skirtVertics.begin(); f != skirtVertics.end(); f++ )
		{
			if(f->pV == p->pV) {	p->pT = f->pT; break;	}
		}
	}
	
	/* 再根据 terrainBorderVertices 各个顶点的pT值 和 faceVertics 的pT值一一对应，替换所有 terrainBorderVertices 里面对应的 pV 和 pVF 值 。*/
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
	{
		for( std::vector<FindValidVertex>::iterator f = faceVertice.begin(); f != faceVertice.end(); f++ )
		{
			if(f->pT == p->pT) {	p->pV = f->pV; p->pVF = f->pVF; break;	}
		}
	}
	
	/* 根据 terrainBorderVertices 各个顶点的pT值 找出四个角点的位置，与conner的下标一致，保存在 isValid里面, 其余都得-1 */
	for( std::vector<FindValidVertex>::iterator f = terrainBorderVertices.begin(); f != terrainBorderVertices.end(); f++ )
	{
		if((f->pT[0] == 0.0) && (f->pT[1] == 0.0))  { f->isValid = 0;}
		if((f->pT[0] == 0.0) && (f->pT[1] == 1.0))  { f->isValid = 1;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 1.0))  { f->isValid = 2;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 0.0))  { f->isValid = 3;}
	}
	
	/* 再还原回 terrainBorderLineLoop */
	terrainBorderLineLoop = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
		terrainBorderLineLoop->push_back(p->pV);

	/*  计算边界上每一个点的经纬度映射。 */
	borderPMs.clear();
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
	{
		ProjectionMap pm;
		getProjectionMapByVertex(&pm, q->pV);
		borderPMs.push_back(pm);	
	}
}

/* 加大裙边的高度 2013-01-04 */
void AddSkirtHeight()
{
	for( std::vector<FindValidVertex>::iterator f= skirtVertics.begin(); f != skirtVertics.end(); f++ )
	{
		osg::Vec3 v = f->pV;
		Geodesy gd;
		Carton  ct;	
		ct.x = v[0]; ct.y = v[1]; ct.z = v[2];	CartToGeod(&ct, &gd); 	
		gd._elevation -= (skirtLength * 0.75);
		GeodToCart(&gd, &ct);
		f->pV.set(ct.x, ct.y, ct.z);
		terrainVertices[f->idx].pV = f->pV;
	}

	/* 改了高度的也要重新写地景 */
	CurrentTerrainIsCut = 1;
}



/* 根据csmList，调整faceVertice各个点的高度，和 coords各个顶点的高度 */
void TuneElevationsByCSMList(osg::ref_ptr<osg::Vec3Array> coords)
{
	if(csmList.size())
	{
		/* 先计算包围球的远近关系*/
		for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = csmList.begin(); c != csmList.end(); c++)
		{
			osg::ref_ptr<CorrectionSpaceMesh> csm = *c;

			/* 计算两个包围球球心之间的距离，和两个包围球半径之和 */	
			float sumRadius = terrainBS.radius() + csm->radius;
			osg::Vec3 Vdis = terrainBS.center() - csm->center;
			float centerDis = Vdis.length();
			if(centerDis < sumRadius) validCSM.push_back(csm);
		}
		if(validCSM.size() == 0) return;
			
#if 1		//计算包围球的办法有点缺陷。如果一个长条的地形就不适合了。
		/* Add by Galen 2013-03-15 采用新办法，计算csmList极点与地形是否相交。*/
		osg::BoundingBox csmListBB;	csmListBB.init();
		for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = csmList.begin(); c != csmList.end(); c++)
		{
			osg::ref_ptr<CorrectionSpaceMesh> csm = *c;
			osg::BoundingBox cBB = csm->getCSMGeode()->getBoundingBox();
			csmListBB.expandBy(cBB);
		}
		
		osg::BoundingBox terBB = pCurrTerrainGeode->getBoundingBox();
		if(terBB.intersects(csmListBB) == false) {	validCSM.clear();	return;	}
		/* Add by Galen 2013-03-15 END */
#endif

		/* 改了高度的也要重新写地景 */
		CurrentTerrainIsCut = 1;
		
		/* 先调整所有面上的点的高度 */
		for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
		{
			for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = validCSM.begin(); c != validCSM.end(); c++)
			{
				osg::ref_ptr<CorrectionSpaceMesh> csm = *c;
				osg::Vec3 fv = f->pV;
				osg::Vec3 Vdis = fv - csm->center;
				float centerDis = Vdis.length();
				if(centerDis < csm->radius)
				{
					csm->correctPoint(f->pV, true);
					terrainVertices[f->idx].pV = f->pV;
				}
			}
		}

		/* 再调整所有裙边上的点的高度 */
		for( std::vector<FindValidVertex>::iterator f= skirtVertics.begin(); f != skirtVertics.end(); f++ )
		{
			for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = validCSM.begin(); c != validCSM.end(); c++)
			{
				osg::ref_ptr<CorrectionSpaceMesh> csm = *c;
				osg::Vec3 fv = f->pV;
				osg::Vec3 Vdis = fv - csm->center;
				float centerDis = Vdis.length();
				if(centerDis < csm->radius)
				{
					csm->correctPoint(f->pV, true);

					osg::Vec3 v = f->pV;
					Geodesy gd;
					Carton  ct;	
					ct.x = v[0]; ct.y = v[1]; ct.z = v[2];	CartToGeod(&ct, &gd); 	gd._elevation -= skirtLength;
					GeodToCart(&gd, &ct);
					f->pV.set(ct.x, ct.y, ct.z);
					terrainVertices[f->idx].pV = f->pV;
				}
			}
		}

		int index = 0;
		for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
		{
			*r = (terrainVertices[index++].pV);
		}
	}
}

/* 在faceVertice集合中找出 角部 和 边部 的一些信息 */
void FindCornerAndBorderInfoInFaceVertices()
{
	/* 在faceVertice集合中查找四个角的顶点信息，从而得出四个边的信息  */
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
	{
		if((f->pT[0] == 1.0) || (f->pT[0] == 0.0) || (f->pT[1] == 1.0) || (f->pT[1] == 0.0))
		{
			FindBorderVertex fbv;	fbv.fvv.pV = f->pV; fbv.fvv.pVF = f->pVF; fbv.fvv.newIdx = f->newIdx; fbv.fvv.idx = f->idx; fbv.fvv.isValid = 0; fbv.onDetailed = 0; fbv.onCulture = 0;
			checkBorderVertices.push_back(fbv);
			f->isValid = 2; // 用来表示面上的点是border上的点
		}else{
			FindValidVertex fvv;	fvv.pV = f->pV; fvv.pVF = f->pVF; fvv.newIdx = f->newIdx; fvv.idx = f->idx; fvv.isValid = 0;
			terrainInnerVertices.push_back(fvv);
		}

		if((f->pT[0] == 0.0) && (f->pT[1] == 0.0))  { conner[0].pV = f->pV; conner[0].pVF = f->pVF; conner[0].pT = f->pT; conner[0].idx = f->idx; conner[0].newIdx = 0; conner[0].isValid = 0;}
		if((f->pT[0] == 0.0) && (f->pT[1] == 1.0))  { conner[1].pV = f->pV; conner[1].pVF = f->pVF; conner[1].pT = f->pT; conner[1].idx = f->idx; conner[1].newIdx = 1; conner[1].isValid = 0;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 1.0))  { conner[2].pV = f->pV; conner[2].pVF = f->pVF; conner[2].pT = f->pT; conner[2].idx = f->idx; conner[2].newIdx = 2; conner[2].isValid = 0;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 0.0))  { conner[3].pV = f->pV; conner[3].pVF = f->pVF; conner[3].pT = f->pT; conner[3].idx = f->idx; conner[3].newIdx = 3; conner[3].isValid = 0;}
	}

	/* 计算将四个角外延的一个角集合 */	
	calcConnerEx();

	/* 计算由四个角组成的简单四条边 */
	border[0].pS = conner[0].pVF; border[0].pE = conner[1].pVF;
	border[1].pS = conner[1].pVF; border[1].pE = conner[2].pVF;
	border[2].pS = conner[2].pVF; border[2].pE = conner[3].pVF;
	border[3].pS = conner[3].pVF; border[3].pE = conner[0].pVF;

	/* 计算一个简单边环，这是在做重构的时候使用的 */
	terrainBoundLoop = new osg::Vec3Array;
	for (unsigned int ic=0;ic<4;ic++) terrainBoundLoop->push_back(conner[ic].pVF);
		terrainBoundLoop->push_back(conner[0].pVF);
}

/* 根据点的信息，确定面的信息。如果组成面的三个点都是面上的点，则这个面就是地景面，否则就是裙边面。 */
void SplitFaceInfoByVertex(osg::Geometry *g)
{	
	for (unsigned int ipr=0; ipr< g->getNumPrimitiveSets(); ipr++) 
	{
		osg::ref_ptr<osg::PrimitiveSet> prset = g->getPrimitiveSet(ipr);
		switch( prset->getType() )
		{
			case osg::PrimitiveSet::DrawArraysPrimitiveType:
			case osg::PrimitiveSet::DrawArrayLengthsPrimitiveType:
			case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
			case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
				break;
			case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
				{
					osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
					if( gdui )
					{
						switch(prset->getMode() )
						{
							case osg::PrimitiveSet::TRIANGLES:
								for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
								{
									unsigned int a,b,c;
									a = gdui->getElement(i);
									b = gdui->getElement(i+1);
									c = gdui->getElement(i+2);				
									if(checkElementUI(skirtVertics, a, b, c))
									{
										faceElement.push_back(a);
										faceElement.push_back(b);
										faceElement.push_back(c);
									}else{
										skirtElement.push_back(a);
										skirtElement.push_back(b);
										skirtElement.push_back(c);
									}
								}
								break;
							default:
								break;
						}
					}
				}
				break;
		}	//case
	}	//for
}

/* 把一个LINE_LOOP，转换成 裙边 的geometry */
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToSkirtFace(osg::ref_ptr<osg::Vec3Array> gvx)
{
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	
	/* 发现个问题，如果线环只有四个点，也就是只有一个三角形的时候，德罗尼三角化会出错。下面分支单独计算。2012-10-10*/
	if(gvx->size() >= 4)
	{
		/* 按输入线环造一个一倍长的线环，*/
		osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;
		for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
		{
			coords->push_back( *pvx );
			osg::Vec3d gI = cTog( *pvx );	gI.z() -= SKIRT_HEIGHT;
			osg::Vec3 cI = gToc( gI );
			coords->push_back( cI );
		}
		gm->setVertexArray(coords);

		/* 再插入面 */
		/* 计算 PrimitiveSet，不同方向要有不同的排列顺序。*/
		osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
		unsigned int left_u = 0, left_d = 1;
		unsigned int right_u, right_d;
		int Direction = 1;
		for(unsigned int ii = 2; ii < coords->size(); ii+=2)
		{
			right_u = ii; right_d = ii + 1;
			if(Direction == 1)
			{
				dui->push_back(left_u); dui->push_back(right_u); dui->push_back(left_d);
				dui->push_back(right_d); dui->push_back(left_d); dui->push_back(right_u);
			}
			if(Direction == 0)
			{
				dui->push_back(left_d); dui->push_back(right_u); dui->push_back(left_u); 
				dui->push_back(right_u); dui->push_back(left_d); dui->push_back(right_d); 
			}
			left_u = right_u; left_d = right_d;
		}
		gm->addPrimitiveSet(dui.get());

	}else{
		gm = NULL;
	}

	return gm.release();
}


/* 把一个LINE_LOOP，转换成 三角面型 的geometry */
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToTranglesFace(osg::ref_ptr<osg::Vec3Array> gvx)
{
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	
	/* 发现个问题，如果线环只有四个点，也就是只有一个三角形的时候，德罗尼三角化会出错。下面分支单独计算。2012-10-10*/
	if(gvx->size() > 4)
	{
		/* 按输入线环造两个相同的线环*/
		osg::ref_ptr<osg::Vec3Array> perimeter1 = new osg::Vec3Array;
		osg::ref_ptr<osg::Vec3Array> perimeter2 = new osg::Vec3Array;
		for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
		{
			perimeter1->push_back( *pvx );
			perimeter2->push_back( *pvx );
		}

		/* 将其中一个线环做德罗尼三角化，另一个线环设为德罗尼约束。 */
		osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(perimeter1.get());
		osg::ref_ptr<ArealConstraint> tdc = new ArealConstraint;
		tdc->setVertexArray( perimeter2.get() );
		tdc->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,perimeter2->size()) );
		dt->addInputConstraint( tdc.get() );
		dt->triangulate();
		dt->removeInternalTriangles( tdc.get() );

		if (tdc.get()->getinteriorTrisSize()>0)
		{
			osg::ref_ptr<osg::Vec3Array> arpts=tdc.get()->getPoints(perimeter1.get());
			gm->setVertexArray(arpts.get());
			gm->addPrimitiveSet(tdc.get()->getTriangles());
		}else
			gm = NULL;
	}else{
		gm->setVertexArray(gvx.get());
		gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES,0,3));
		osg::Vec3Array::iterator r = gvx->begin();
		osg::Vec3 x, y, z;	x = *r++;	y = *r++;	z = *r++;
		double area = getArea(x, y, z);
		if(area < 10.0)	gm = NULL;
	}

	return gm.release();
}

/* 从一个dclist里面找出所有有效点（就是参与组成线环的点）的集合，*/
osg::ref_ptr<osg::Vec3Array> FindValidVertexInDC(osg::ref_ptr<ArealConstraint> dc)
{
	osg::ref_ptr<osg::Vec3Array> vtdc = dynamic_cast<osg::Vec3Array *> (dc->getVertexArray());
	osg::ref_ptr<osg::Vec3Array> vtx = new osg::Vec3Array;
	if(dc->getNumPrimitiveSets() == 0) return NULL;
	osg::ref_ptr<osg::PrimitiveSet> prset = dc->getPrimitiveSet(0);

	switch( prset->getType() )
	{
		case osg::PrimitiveSet::DrawArraysPrimitiveType:
			{
				osg::DrawArrays* gdui = dynamic_cast<osg::DrawArrays*>(prset.get());
				if( gdui )
				{
					switch(prset->getMode() )
					{
						case osg::PrimitiveSet::LINE_LOOP:
							for(unsigned int i=0; i<gdui->getNumIndices(); i++)
							{
								unsigned int a = gdui->index(i);
								vtx->push_back(*(vtdc->begin()+a));
							}
							break;
						default:
							break;
					}
				}
			}
			break;
		case osg::PrimitiveSet::DrawArrayLengthsPrimitiveType:
			break;
		case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
			{
				osg::DrawElementsUShort *gdui = dynamic_cast<osg::DrawElementsUShort *>(prset.get());
				if( gdui )
				{
					switch(prset->getMode() )
					{
						case osg::PrimitiveSet::LINE_LOOP:
							for(unsigned int i=0; i<gdui->getNumIndices(); i++)
							{
								unsigned int a = gdui->getElement(i);
								vtx->push_back(*(vtdc->begin()+a));
							}
							break;
						default:
							break;
					}
				}
			}
			break;
		case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
			{
				osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
				if( gdui )
				{
					switch(prset->getMode() )
					{
						case osg::PrimitiveSet::LINE_LOOP:
							for(unsigned int i=0; i<gdui->getNumIndices(); i++)
							{
								unsigned int a = gdui->getElement(i);
								vtx->push_back(*(vtdc->begin()+a));
							}
							break;
						default:
							break;
					}
				}
			}
			break;
		case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
			{
				osg::DrawElementsUByte *gdui = dynamic_cast<osg::DrawElementsUByte *>(prset.get());
				if( gdui )
				{
					switch(prset->getMode() )
					{
						case osg::PrimitiveSet::LINE_LOOP:
							for(unsigned int i=0; i<gdui->getNumIndices(); i++)
							{
								unsigned int a = gdui->getElement(i);
								vtx->push_back(*(vtdc->begin()+a));
							}
							break;
						default:
							break;
					}
				}
			}
			break;
	}	//case
	
	return vtx.get();
}


/* 把一个dc，转换成 三角面型 的geometry */
osg::ref_ptr<osg::Geometry> ConvertFromDCToTranglesFace(osg::ref_ptr<ArealConstraint> dc)
{
	osg::ref_ptr<osg::Vec3Array> vtx = FindValidVertexInDC(dc);
	if( vtx == NULL ) return NULL;

	/* 使vtx线环保证首尾相连接 */
	if (!vtx->empty() && vtx->front()!=vtx->back())
		vtx->push_back(vtx->front());

	osg::ref_ptr<osg::Geometry> gm = ConvertFromLINELOOPToTranglesFace(vtx);
	return gm.release();
}

#ifdef ADD_SKIRT
	/* 创建裙边 */
void CreateSkirt(osg::ref_ptr<osg::Vec3Array> validCoords, osg::ref_ptr<osg::DrawElementsUInt> dui, int mode)
{
	/* 首先初始化全局变量 */
	cultSkirtVertice.clear();
	cultSkirtDrawElements.clear();
	
	/* 分与地形边界有交点 和 与地形边界无交点 两种情况，分别计算裙边 */
	if (intersectlist.empty())
	{
		if(mode == 1)
		{
			//将原来裙边的点也加入到点集合当中，并确定裙边点新序号newIdx
			int startIdx = validCoords->size();
			for( std::vector<FindValidVertex>::iterator p = skirtVertics.begin(); p != skirtVertics.end(); p++ )
			{
				validCoords->push_back(p->pV); 
				p->newIdx = startIdx++;
			}
	
			//原来做德罗尼之前的面上的点集合，都要在有效点集合validVertices里面查找出在做了德罗尼之后的新的idx;
			for( std::vector<FindValidVertex>::iterator p = faceVertice.begin(); p != faceVertice.end(); p++ )
			{
				for( std::vector<FindValidVertex>::iterator pp = validVertices.begin(); pp != validVertices.end(); pp++ )
	  			{
	  				if(p->pV == pp->pV)
	  				{
	  					p->newIdx = pp->newIdx; break;
	  				}
	  			}
			}
	 
			//将裙边面元素集skirtElement里面老的idx替换成新的newIdx信息
			for( std::vector<unsigned int>::iterator p = skirtElement.begin(); p != skirtElement.end(); p++ )
			{
				int flag = 0;
				for( std::vector<FindValidVertex>::iterator ppp = faceVertice.begin(); ppp != faceVertice.end(); ppp++ )
				{
					if (*p==ppp->idx) 
					{
						*p = ppp->newIdx;
						flag = 1;
						break;
					}
				}
	
				if(flag) continue;
				for( std::vector<FindValidVertex>::iterator pp = skirtVertics.begin(); pp != skirtVertics.end(); pp++ )
				{
					if (*p==pp->idx) 
					{
						*p = pp->newIdx;
						break;
					}
				}
			}
	
			//再增加所有裙边skirt的面
			for( std::vector<unsigned int>::iterator p = skirtElement.begin(); p != skirtElement.end(); p++ )
				dui->push_back(*p);
				
			/* 先把裙边的点也加入到validVertices集合里面 */
			for( std::vector<FindValidVertex>::iterator p = skirtVertics.begin(); p != skirtVertics.end(); p++ )
			{
				FindValidVertex fvv;
				fvv.pV = p->pV; fvv.pVF = p->pVF; fvv.pT = p->pT; fvv.idx = p->idx; fvv.newIdx = p->newIdx; fvv.isValid = p->isValid;
				validVertices.push_back(fvv);
			}
	
			/* 根据边界段端点的有效性，来决定裙边是不是可以画出来。*/
			//找出所有的有效面。存放在 validElement 里面
			osg::ref_ptr<osg::DrawElementsUInt> validdui = dui;
			validElement.clear();
			for(unsigned int i=0; i<validdui->getNumIndices(); i+=3)
			{
				unsigned int a,b,c;
				a = validdui->getElement(i);
				b = validdui->getElement(i+1);
				c = validdui->getElement(i+2);
				if(checkTriOnSkirt(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF))
				{
					validElement.push_back(a);
					validElement.push_back(b);
					validElement.push_back(c);
				}
			}
			//把只包含合法三角形的validElement里面的元素加进来
			dui->clear();
			for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
				dui->push_back(*p);
		}
	}
	else
	{
		osg::ref_ptr<osg::Geometry > dtGeometry = new osg::Geometry;    
		dtGeometry->setVertexArray( validCoords.get() );
		dtGeometry->addPrimitiveSet( dui );
		TerrainProcessingUtils::computePerimeterAll( dtGeometry.get(), perimeterList);
	}

	/* 制作裙边之前，检查一下线环的面积是否合适。*/
	{
		std::vector<osg::ref_ptr<osg::Vec3Array>> temp;
		for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perimeterList.begin(); q != perimeterList.end(); q++ )
			temp.push_back( *q );
		perimeterList.clear();
		for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= temp.begin(); q != temp.end(); q++ )
		{
			double area = GetAreaOfLineloop(*q);
			if(area > AREA_LIMIT)	perimeterList.push_back(*q);
		}
	}
	if(perimeterList.size() == 0)	return;

	/* 制作裙边 */
	makeSkirt(perimeterList, validCoords->size());

	/* 由于在makeSkirt()里面有可能修改了validVertices里面的值，所以需要重新创建validCoords */
	validCoords->clear();
	for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
		validCoords->push_back(pp->pV);

	/* 最后加入新做好的裙边点 */
	for( std::vector<FindValidVertex>::iterator p = cultSkirtVertice.begin(); p != cultSkirtVertice.end(); p++ )
		validCoords->push_back(p->pV);
}

#endif	//ADD_SKIRT


/* 处理每一个Geometry */
osg::ref_ptr<osg::Geometry> processGeometry( osg::Geometry *g, const GeodeList &features )
{
	elevationOfLake = -9999.99;
	pCurTerrain = g;
	pCurrTerrainGeode = new osg::Geode;
	pCurrTerrainGeode->addDrawable(g);
	osg::ref_ptr<osg::Vec3Array> coords    = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	osg::ref_ptr<osg::Vec2Array> texCoords = dynamic_cast<osg::Vec2Array *>(g->getTexCoordArray(0));
	osg::ref_ptr<osg::Vec3Array> norCoords = dynamic_cast<osg::Vec3Array *>(g->getNormalArray() );

	/* 地景面有效性检查 */
	if( !coords.valid() )    return NULL;
	if( coords->size() == 0) return NULL;
    
	/* 计算当前地景块的包围球 */
	terrainBS = g->getBound();
    
	/* 第一步，将裙边skirt点面从原Geometry中排除出去。 */
	/* 新创建一个Geometry，不包含输入参数g里面的裙边skirt点和面。*/
	/* 初始化所有用到的数组 */
	InitVectors();

	/* 读入所有顶点信息和纹理信息到terrainVertices数组, 并计算出投影到海平面的各个点的值。 */
	int j = 0;
	osg::Vec2Array::const_iterator pp = texCoords->begin();
	osg::Vec3Array::const_iterator pn = norCoords->begin();
	for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end() && pp!=texCoords->end() && pn!=norCoords->end(); r++ )
	{
		FindValidVertex fvv;
		fvv.pV = *r; fvv.pT = *pp; fvv.pN = *pn; fvv.idx = j; fvv.newIdx = 0; j++; pp++; pn++;
		terrainVertices.push_back(fvv);
	}
	getFlatVertex(&terrainVertices);
    
	/* 将整个地景面上的所有点分为两个部分，一部分是地景面上的点，一部分是裙边上的点。*/
	SplitTerrainVertexToFaceAndSkirt();

	/* 如果没有裙边，则不是地景，有可能是 features 文件中不被识别的gnode */
	if (skirtVertics.empty()) return NULL;

#ifdef ADD_SKIRT
	/* 计算现有裙边的高度 */
	Geodesy fgeo, sgeo;
	Carton  fcat, scat;
	fcat.x = faceVertice[0].pV[0];	fcat.y = faceVertice[0].pV[1];	fcat.z = faceVertice[0].pV[2];
	scat.x = skirtVertics[0].pV[0];	scat.y = skirtVertics[0].pV[1];	scat.z = skirtVertics[0].pV[2];
	assert(skirtVertics[0].pT==faceVertice[0].pT);
	CartToGeod(&fcat, &fgeo);	CartToGeod(&scat, &sgeo);
	skirtLength = fgeo._elevation - sgeo._elevation;
#endif

	/* 调整裙边的点，加大裙边的高度 2013-01-04 */
#ifdef	MAKE_HIGH_SKIRT
	skirtLength *= 4;
	AddSkirtHeight();
#endif

	/* 根据csmList，调整faceVertice和coords各个顶点的高度 */
	TuneElevationsByCSMList(coords);


	/* 初始化纹理坐标保存器。 */
	osgTDS::TexCoordPreserver tcp(&terrainVertices);
    	
	/* 在faceVertice集合中找出 角部 和 边部 的一些信息 */
	FindCornerAndBorderInfoInFaceVertices();

	/* 将边界线环上的所有点，由裙边点，都替换成面上的点 */
	ChangeBorderVertexFromSkirtToFace(g);
	
	/* 根据点的信息，确定面的信息。如果组成面的三个点都是面上的点，则这个面就是地景面，否则就是裙边面。 */
	SplitFaceInfoByVertex(g);
	
	// 第三步，做德罗尼，只用面上的点，排除裙边的点。
	// 先创建一个Geode，只包含面上的三角形，点集合就全包括，而且这些点都转换成投影点，投影到海平面。(也就是点转换成经纬高坐标系后，高度值都改为0.0)，用在做德罗尼。
	osg::ref_ptr<osg::DrawElementsUInt> facedui=new osg::DrawElementsUInt(GL_TRIANGLES);
	for(std::vector<unsigned int>::iterator p = faceElement.begin(); p!=faceElement.end(); p++)
		facedui->push_back(*p);
	osg::ref_ptr<osg::Geometry > faceGeometry = new osg::Geometry;
	osg::ref_ptr<osg::Vec3Array> coords_Flat = getCoordsFlat(terrainVertices);
	faceGeometry->setVertexArray( coords_Flat.get() );
	faceGeometry->addPrimitiveSet( facedui.get() );
	osg::ref_ptr<osg::Geode> flatTerrainGd = new osg::Geode;
	flatTerrainGd->addDrawable(faceGeometry.get());


	IntersectCt = 0;
	std::vector<osg::ref_ptr<osg::Geometry>> gmList;
	std::vector< osg::ref_ptr<osg::Vec3Array> > roadList;
	osg::ref_ptr<osg::Geode> lakeGd = new osg::Geode;

	int ic;
	f_idx = 0;
	bool hasDetailed = false;
	bool hasLakes = false;
	bool isAirport = false;
	for( GeodeList::const_iterator p = features.begin(); p != features.end(); p++ )
	{
		subIntersectCt = 0;
		gmList.clear();
		featurePMs.clear();

		/* 先把地形方块四个角清零，预备计算当前features与地形方块四个角的相交情况。*/
		for (ic=0;ic<4;ic++)	conner[ic].isValid = 0;

		osg::Geode* f = (*p)._geode.get();
		
		/* 下面定义的flatFeatureGd, 是文化信息上所有点都做了投影的Geode。*/
		osg::ref_ptr<osg::Geode> flatFeatureGd = new osg::Geode;
		for( unsigned int ii = 0; ii < f->getNumDrawables(); ii++ )
		{
			osg::ref_ptr<osg::Geometry> g = dynamic_cast<osg::Geometry *>(f->getDrawable(ii));
			osg::ref_ptr<osg::Geometry> fg = new osg::Geometry(*g, osg::CopyOp::DEEP_COPY_ALL);
			std::string fname;
			if (g->getName().empty()==false)  fname = g->getName(); else fname = f->getName();

			
			/* 如果当前名称是Lake，那么获取湖面的高度 */
			if(fname.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			{
				if(elevationOfLake < 0.0)
				{
					osg::ref_ptr<osg::Vec3Array> ver =  dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
					osg::Vec3Array::iterator q = ver->begin();
					
					/* 取其中一个点的高度即可。*/
					Geodesy geo;
					Carton  cat;
	
					cat.x = q->x(); cat.y = q->y(); cat.z = q->z();	CartToGeod(&cat, &geo);	
					double currLakeElev = geo._elevation;
					
					if(currLakeElev > elevationOfLake) elevationOfLake = currLakeElev;
				}
			}

			/* 如果当前名称是Airport，那么获取机场的高度 */
			isAirport = false;
			m_elevOfAirport = 0.0;
			if((fname.substr(0, strlen(DEFINE_AIRPORT_STRING)) == DEFINE_AIRPORT_STRING) || (fname.substr(0, 7) == "AIRPORT"))
			{
				SplitName(fname);
				if(NameParts.size() >= 2)	m_elevOfAirport = atof(NameParts[1].c_str());
				isAirport = true;
			}

			/* 获取顶点个数。*/
			osg::Vec3Array *vert =  dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
			if(vert == NULL)	continue;
			int vert_size =  vert->size();

			/* 假如Geometry类型为LINE_LOOP，要把LINE_LOOP恢复成面的类型。 */
			osg::ref_ptr<osg::PrimitiveSet> ps = g->getPrimitiveSet(0);
			if(ps->getMode() == osg::PrimitiveSet::LINE_LOOP)
			{
				osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
				fg = ConvertFromLINELOOPToTranglesFace(gvx);
				if(fg == NULL)	continue;
			}
			
			/*把fg所有顶点都做投影化处理，投影到海平面。*/
			makeGeometryFlat(fg);
			flatFeatureGd->addDrawable(fg);

			/* 重构features之前，还是要先算一下地景四个角与features的相交情况，
			   因为：出现了地景一个角，与features相交，但却与重构后的features不相交 的情况。造成错误。
			   地景方块的四个角点与features是否有交点，有交点就让isVaild = 1  */
			osg::Vec3 isect;
			for (ic=0;ic<4;ic++)
			{
				if( TerrainProcessingUtils::getIntersect( flatFeatureGd, conner[ic].pVF, isect ))	conner[ic].isValid += 1; // 表示角点跟 feature 相交
			}

			/* 查看面上所有的点是否跟 feature 相交。 */
			int TotalIntersectVerts = 0;
			for( std::vector<FindValidVertex>::iterator q= faceVertice.begin(); q != faceVertice.end(); q++ )
			{
				if (TerrainProcessingUtils::getIntersect( flatFeatureGd, q->pVF, isect ))	
				{	
					subIntersectCt++;	
					TotalIntersectVerts++;	
					
					/* 如果文化信息面是机场面，*/
					if((isAirport == true) && (m_elevOfAirport != 0.0))
					{
						Geodesy geo;
						Carton  cat;
					
						cat.x = q->pVF.x(); cat.y = q->pVF.y(); cat.z = q->pVF.z();
						CartToGeod(&cat, &geo);	geo._elevation = m_elevOfAirport; GeodToCart(&geo, &cat);
						q->pV.x() = cat.x; q->pV.y() = cat.y; q->pV.z() = cat.z;
						
						/* 设置更新高程的标志变量 */
						UpdateTheElevationInTerrain = 1;
					}
				}
			}
			
#ifdef ADD_SKIRT
			/* 补充一下，再查看裙边上所有的点是否跟 feature 相交。*/
			for( std::vector<FindValidVertex>::iterator q= skirtVertics.begin(); q != skirtVertics.end(); q++ )
			{
				if (TerrainProcessingUtils::getIntersect( flatFeatureGd, q->pVF, isect ))	
				{	
					/* 如果文化信息面是机场面，*/
					if((isAirport == true) && (m_elevOfAirport != 0.0))
					{
						Geodesy geo;
						Carton  cat;
					
						cat.x = q->pVF.x(); cat.y = q->pVF.y(); cat.z = q->pVF.z();
						CartToGeod(&cat, &geo);	geo._elevation = m_elevOfAirport - skirtLength; GeodToCart(&geo, &cat);
						q->pV.x() = cat.x; q->pV.y() = cat.y; q->pV.z() = cat.z;

						/* 设置更新高程的标志变量 */
						UpdateTheElevationInTerrain = 1;
					}
				}
			}
#endif
			
			/* 如果地景面上所有的点都和文化信息面相交，那么可以把文化信息面简单处理一下。直接覆盖这个地景面就可以。*/
			osg::ref_ptr<osg::Geometry> gg;
			if(TotalIntersectVerts == faceVertice.size())
			{
				gg = new osg::Geometry;
				osg::ref_ptr<osg::Vec3Array> vtx = new osg::Vec3Array;
				vtx->push_back(connerEx[0].pVF);
				vtx->push_back(connerEx[1].pVF);
				vtx->push_back(connerEx[2].pVF);
				vtx->push_back(connerEx[3].pVF);
				vtx->push_back(connerEx[0].pVF);
				vtx->push_back(connerEx[2].pVF);
				gg->setVertexArray(vtx);
				gg->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES,0,6));
			}else
				gg = new osg::Geometry(*fg, osg::CopyOp::DEEP_COPY_ALL);

			/* Add by Galen 2013-02-21 */
			/* 如果当前文化信息是湖泊，那么要对临近湖边的顶点做一些收集，以便后面做平滑处理。*/
			if(fname.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			{
				getVextexListNearToCulture(flatTerrainGd.get(), gg.get());
			}
			/* Add by Galen 2013-02-21 End */

			/* 重构features, 把 features 上不跟地景相交的面删掉，处理路跟地景有多次相交的情况 */
			/* 这里用投影后的地景面 和 投影后的文化信息features面（features的投影在函数内做。）做重构，精确度能更高一些。*/
			if (reConstructFeature(flatTerrainGd.get(), gg.get(), terrainBoundLoop.get()))
			{
				std::string name;
				if (g->getName().empty()==false)  name = g->getName(); else name = f->getName();
				gg->setName(name);
				
				gmList.push_back(gg);
			}
		}	//geode
		
		/* 对于每个geode的所有geometry, 求出外围环点集合。*/
		for( std::vector<osg::ref_ptr<osg::Geometry> >::iterator qm= gmList.begin(); qm != gmList.end(); qm++ )
		{
			/* 对所有边界的点进行清理，如果边界点和细节纹理区块或者岛相交，则这个边界点是有效边界点。*/
			std::string gmName = qm->get()->getName();
			osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
			gnode->addDrawable(qm->get());
			int nDetail, nCulture; nDetail = nCulture = 0; 
			if (isDetailedtxt(gmName))
			{
				for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
				{
					osg::Vec3 isect;
					if (TerrainProcessingUtils::getIntersect( gnode, p->fvv.pVF, isect ))	
					{	
						FindValidVertex fvv;	fvv.pV = p->fvv.pV; fvv.pVF = p->fvv.pVF; fvv.newIdx = p->fvv.newIdx; fvv.idx = p->fvv.idx; fvv.isValid = 0; validBorderVertices.push_back(fvv);
						nDetail++; p->onDetailed = 1;
					}
				}
			}
			else{
				for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
				{
					osg::Vec3 isect;
					if (TerrainProcessingUtils::getIntersect( gnode, p->fvv.pVF, isect ))	
					{	
						nCulture++;	p->onCulture = 1;
					}
				}
			}			

			/* 计算线环的集合 */
			std::vector<osg::ref_ptr<osg::Vec3Array> > perList;
			if (qm->get()->getPrimitiveSet(0)->getMode() == osg::PrimitiveSet::LINE_LOOP)
			{
				osg::ref_ptr<osg::Vec3Array> perimeter = new osg::Vec3Array;
				osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(qm->get()->getVertexArray());
				for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
					perimeter->push_back( *pvx );
	
				perList.push_back(perimeter);
			}
			else
			{
				TerrainProcessingUtils::computePerimeterAll( qm->get(), perList);
			}

			for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perList.begin(); q != perList.end(); q++ )
			{
				/* 先把当前perlist线环都三角化成一个面，为了跟边界段集合做相交检查 */				
				osg::ref_ptr<osg::Geometry> fg = ConvertFromLINELOOPToTranglesFace(*q);
				if(fg == NULL)	continue;
				if((qm->get()->getName().substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING))
				{
					lakeGd->addDrawable(fg);
				}

				/* 查看 terrainBorderVertices 的所有边界点是否与 当前perlist线环面 相交，如果相交，那么 newIdx 项 = 1，否则 = 0 */
				borderIntersectLS.clear();
				borderHalfIntersectLS.clear();
				osg::ref_ptr<osg::Geode> currPerlistGd = new osg::Geode;
				currPerlistGd->addDrawable(fg.get());
				int currentBorderIdx = -1;
				for( std::vector<FindValidVertex>::iterator bq= terrainBorderVertices.begin(); bq != terrainBorderVertices.end(); bq++ )
				{
					osg::Vec3 isect;
					if (TerrainProcessingUtils::getIntersect( currPerlistGd, bq->pVF, isect ))	{	bq->newIdx = 1; } else {	bq->newIdx = 0; }
					if(bq != terrainBorderVertices.begin())
					{
						std::vector<FindValidVertex>::iterator bp = bq - 1;
						if(bp->isValid == 0) currentBorderIdx = 3;
						if(bp->isValid == 1) currentBorderIdx = 0;
						if(bp->isValid == 2) currentBorderIdx = 1;
						if(bp->isValid == 3) currentBorderIdx = 2;
						
						LineSeg ls;
						ls.pS = bp->pVF;   ls.idxS = bp->idx;	ls.validS = bp->newIdx;
						ls.pE = bq->pVF;   ls.idxE = bq->idx;	ls.validE = bq->newIdx;	ls.border = currentBorderIdx;
						borderIntersectLS.push_back(ls);
	
						/* 下面一个判断仅供设断点用 */
						if(((bp->newIdx == 1)&&(bq->newIdx == 0))||((bp->newIdx == 0)&&(bq->newIdx == 1)))
						{
							borderHalfIntersectLS.push_back(ls);
						}
					}
				}

				for( std::vector<LineSeg>::iterator bi= borderIntersectLS.begin(); bi != borderIntersectLS.end(); bi++ )
				{
					if(bi->border == -1) bi->border = currentBorderIdx; else break;
				}
				
				/* 查看 terrainInnerVertices 的所有内部点是否与 当前perlist线环面 相交，如果相交，那么 AllInnerVertexIntersection 项 = true，否则 = false */
				AllInnerVertexIntersection = false;
				AllInnerVertexNoIntersection = false;
				int TotalInnerIntersectVerts = 0;
				osg::Vec3 isect;
				bool isDetail = isDetailedtxt(gmName);
				for( std::vector<FindValidVertex>::iterator qq= terrainInnerVertices.begin(); qq != terrainInnerVertices.end(); qq++ )
				{
					if (TerrainProcessingUtils::getIntersect( currPerlistGd, qq->pVF, isect ))	
					{	
						TotalInnerIntersectVerts++;	
						
						/* 如果内部点与perlist线环面相交，那么这个内部点可以作为特殊点存起来。*/
						if(isDetail == false)	innerVerticesInterWithCulture.push_back(qq->pVF);
					}
				}
				if(TotalInnerIntersectVerts == terrainInnerVertices.size()) AllInnerVertexIntersection = true;
				if(TotalInnerIntersectVerts == 0) AllInnerVertexNoIntersection = true;
	
				/* 增加几何信息与地景边界的交叉点 */
				needDoSkirt = true;
				osg::ref_ptr<osg::Vec3Array> pmForCulture = new osg::Vec3Array;
				addIntersectVertex(flatTerrainGd, currPerlistGd, *q, pmForCulture, gmName);
	
				for(std::vector< osg::ref_ptr<osg::Vec3Array> >::iterator mq = multiDc.begin(); mq != multiDc.end(); mq++)
				{
					if (mq->get()->size()>=3) 
					{
						// dc 用来生成细节纹理的面
						std::string name = qm->get()->getName();

						/* 把路, 湖先存起来 */
						if(name.substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)	roadList.push_back(*mq);
						if(name.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)	hasLakes = true;

						osg::ref_ptr<ArealConstraint> dc = new ArealConstraint;
						if (isDetailedtxt(name))
						{
							hasDetailed = true;
						}
						else
						{
							/* 这个外围环集合，是专门用来做裙边用的。*/
							if(needDoSkirt)
							{
								/* 如果有动态草，则不需要做裙边。*/
								if(name.substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) != DEFINE_DYNAMICGRASS_STRING)
								{
									perimeterList.push_back(*mq);
								}
							}
						}

						/* dc 是用来生成文化信息的面 */
						dc.get()->setName(name);
						dc->setVertexArray( mq->get() );
						dc->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0, mq->get()->size()) );
						dclist.push_back( dc.get() );

					}
				}
			}	//PerList
		}	//gmList
		f_idx++;
		
		IntersectCt += subIntersectCt;
	}	//features

	/* 把所有的地形面的顶点，按照湖面调整一下高度 */
	TuneElevationOfLakeSurface(lakeGd, coords);

	/* 完善边界点的集合 */
	/* 在前面循环里面，已经添加所有的“文化信息与地形边界交点”进入了allBorderVertices。*/
	/* 在这一步，把所有地形本身的边界点，也加入到 allBorderVertices */
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
		allBorderVertices.push_back(*q);

	/* 先把所有的地形面的顶点，都集中到totalVertices集合里面 */
	getTotalFaceVertex();

	/* 把路保存到细节纹理区 */
	if(roadList.size() > 0)
	{
		for(std::vector< osg::ref_ptr<osg::Vec3Array> >::iterator q = roadList.begin(); q != roadList.end(); q++)
		{
			/*   第一步    */
			/* 把当前perlist线环都三角化成一个面 */				
			osg::ref_ptr<osg::Geometry> gm = ConvertFromLINELOOPToTranglesFace(*q);
			if(gm == NULL)	continue;

			osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
			osg::ref_ptr<osg::StateSet> st = new osg::StateSet(*(pCurTerrain->getStateSet()));
			osg::Texture2D* t2d = dynamic_cast<osg::Texture2D*>(st->getTextureAttribute(0,osg::StateAttribute::TEXTURE));
			st->setTextureAttributeAndModes(0,t2d,osg::StateAttribute::ON);
			gnode->setStateSet(st.get());
			gnode->setName("RoadDay");
			gnode->addDrawable(gm);
			detailed_gdlist.push_back(gnode); 

			/* 增加纹理坐标 */
			osg::ref_ptr<osg::Vec3Array> validCoords = dynamic_cast<osg::Vec3Array *>(gm->getVertexArray());

			/* 把平面点映射回空间点 */
			for( osg::Vec3Array::iterator r = validCoords->begin(); r != validCoords->end(); r++ )
			{
				osg::Vec3 pSpace = getSpaceVertex(*r);
				*r = pSpace;
			}
			
			osg::ref_ptr<osg::Vec2Array> getc = tcp.getTexCoords( *(validCoords.get()) );
			gm->setTexCoordArray( 0, getc.get() );
			// Restore normal 
			osg::ref_ptr<osg::Vec3Array> nc = tcp.getNormals(  *(validCoords.get()) );
			gm->setNormalArray( nc.get() );
			gm->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
			// Add an overall white color
			osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
			color->push_back(osg::Vec4(1,1,1,1));
			gm->setColorArray(color.get());
			gm->setColorBinding( osg::Geometry::BIND_OVERALL );
				
			/*   第二步    */
			/* 把当前perlist线环制作成一个裙边 */
			osg::ref_ptr<osg::Geometry> gmSkirt = ConvertFromLINELOOPToSkirtFace(*q);
			if(gmSkirt == NULL) continue;

			/* 增加纹理坐标 */
			osg::ref_ptr<osg::Vec3Array> skirtCoords = dynamic_cast<osg::Vec3Array *>(gmSkirt->getVertexArray());

			/* 把平面点映射回空间点 */
			for( osg::Vec3Array::iterator r = skirtCoords->begin(); r != skirtCoords->end(); r++ )
			{
				osg::Vec3 pSpace = getSpaceVertexEx(*r);
				*r = pSpace;
			}

			osg::ref_ptr<osg::Vec2Array> gets = tcp.getTexCoords( *(skirtCoords.get()) );
			gmSkirt->setTexCoordArray( 0, gets.get() );
			// Restore normal 
			osg::ref_ptr<osg::Vec3Array> ns = tcp.getNormals(  *(skirtCoords.get()) );
			gmSkirt->setNormalArray( ns.get() );
			gmSkirt->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
			// Add an overall white color
			osg::ref_ptr<osg::Vec4Array> colorS = new osg::Vec4Array;
			colorS->push_back(osg::Vec4(1,1,1,1));
			gmSkirt->setColorArray(colorS.get());
			gmSkirt->setColorBinding( osg::Geometry::BIND_OVERALL );
			gnode->addDrawable(gmSkirt);

		}
	}


	/* 把细节纹理的dc都单独拿出来 */
	MoveOutDetailTxtDc();

	/* 整理边界点，把纯粹被文化信息覆盖的点加入到无效点集合invalidVertices里面 */
	int numInvalidBorderVertex = 0;
	bool allBorderInvalid = false;

#if 0	//如果德罗尼运算正确的话，这段程序实际上也可以不要。
	for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
	{
		if((p->onDetailed == 0) && (p->onCulture == 1)) 
		{
			invalidVertices.push_back(p->fvv.pVF);
			numInvalidBorderVertex++;
		}
	}
#endif

	if(numInvalidBorderVertex == checkBorderVertices.size()) allBorderInvalid = true;

	/* 把同名的多个dc合并到一起，组合成一个dc */
	/* 把所有具有相同名称的dc合并成为一个，并计算覆盖情况 */
	MergeDclistStepOne(&dclist);

	/* 如果dclist里面只有两个dc, 且两个dc一模一样。那么就舍弃其中一个。*/
	CheckDclistAndMerge(&dclist);

	//如果没有有效的约束生成，及早退出
	if (dclist.empty() && (dcDetail.empty() && dcLakes.empty() ))
	{
		/* 没有有效德罗尼约束生成的话，...... */
		if (IntersectCt < faceVertice.size())			// < 3
		{
			/* 当前返回这块虽然没有被切割, 看看是否被调整了高度. */
			if(validCSM.size())
			{
				/* 还要根据修改后的faceVertice更正原有coords里面对应的顶点值，否则会出现断层 */
				/* 因为存有顶点序号，所以直接按照序号指示修改即可 */
				for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
				{
					int idx = f->idx;	terrainVertices[idx].pV = f->pV;
				}
				int index = 0;
				for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
				{
					*r = (terrainVertices[index++].pV);
				}

				/* 改了高度的也要重新写地景 */
				CurrentTerrainIsCut = 1;
			}
			
#ifdef	MAKE_HIGH_SKIRT
			/* 由于裙边有被加长 更改一下裙边的顶点。2013-01-04 */
			if(1)
			{
				/* 因为存有顶点序号，所以直接按照序号指示修改即可 */
				for( std::vector<FindValidVertex>::iterator f= skirtVertics.begin(); f != skirtVertics.end(); f++ )
				{
					int idx = f->idx;	terrainVertices[idx].pV = f->pV;
				}
				int index = 0;
				for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
				{
					*r = (terrainVertices[index++].pV);
				}
			}
#endif

			return g;
		}
		else 
		{
			CurrentTerrainIsCut = 1;
			return NULL;
		}
	}

	/* 如果地景有变化，设置全局通讯变量。*/
	CurrentTerrainIsCut = 1;

	//如果少于1个面上的点跟feature不相交，不要这块地景
	int VertexNoIntersect = abs((int)faceVertice.size() - IntersectCt );
	if(( VertexNoIntersect < 1) && (!hasDetailed) && (dclist.size() == 1))
		return NULL;

	/***********************************************************************/
	/********  制作细节纹理面  *********************************************/
	if((dcDetail.size() > 0) && (dclist.size() == 0))
	{
		/* 生成纹理信息的面 */
		CreateFaceOfDetailTxtWhenNoCulture(&tcp);

		/* 这里描述了一个状况：文化信息和细节纹理占满了整个地形，按道理应该在前面退出；
		   但是要求出细节纹理之后才能退出，否则地面就缺少了一块细节纹理面了。
		 */ 
		if((VertexNoIntersect < 1 ) && (hasDetailed))
			return NULL;
		
		/* 将dcDetail赋值给dclist */
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
		{
			dclist.push_back(*p);
		} 
		dcDetail.clear();
	}

	/***********************************************************************/
	/********  制作湖面  *********************************************/
	if( (dcLakes.size() > 0) )
	{
		/* 如果满足下面三个条件，说明整个地块都是湖面，要提前退出。*/
		if(( VertexNoIntersect < 1) && (hasLakes) )
		{
			osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
			osg::ref_ptr<osg::StateSet> st = new osg::StateSet(*(pCurTerrain->getStateSet()));
			osg::Texture2D* t2d = dynamic_cast<osg::Texture2D*>(st->getTextureAttribute(0,osg::StateAttribute::TEXTURE));
			st->setTextureAttributeAndModes(0,t2d,osg::StateAttribute::ON);
			gnode->setStateSet(st.get());
			gnode->setName("Lake");
			gnode->addDrawable(g);
			detailed_gdlist.push_back(gnode);
			dcLakes.clear();
			return NULL;
		}
		
		/* 把Lake和Island的dc合并成为一个，并计算覆盖情况 */
		MergeLakeAndIsland(&dclist);

		/* 生成湖面 */
		CreateFaceOfLake(&tcp);
		
		/* 如果有岛的话，把岛从dclist里面去掉 */
		std::vector<osg::ref_ptr<ArealConstraint> > noIslands;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
		{
			std::string name= p->get()->getName();
			if (name.substr(0, strlen(DEFINE_ISLAND_STRING)) != DEFINE_ISLAND_STRING)
				noIslands.push_back(*p);
		}
		dclist.clear();
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= noIslands.begin(); p != noIslands.end(); p++ )
			dclist.push_back(*p);

		/* 将dcLakes 赋值给 dclist */
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
		{
			dclist.push_back(*p);
		} 
		dcLakes.clear();
	}


	/***********************************************************************/
	/********  制作动态草的面  *********************************************/
	if((dclist.size() != 0) )
	{

		/* 检查dclist里面是否有DynamicGrass，如果没有就退出。*/
		int iFlag = 0;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
		{
			if(p->get()->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) == DEFINE_DYNAMICGRASS_STRING)
			{	iFlag = 1; break;	}
		}

		if(iFlag)
		{
			if(CurrentTerrainLodLevel == 4)
			{
				/* 生成动态草的面 */
				CreateFaceOfDynamicGrass(&tcp);
			}
	
			/* 将动态草dc从dclist中去除 */
			std::vector<osg::ref_ptr<ArealConstraint> > tmp_dclist;
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
			{
				if (p->get()->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) != DEFINE_DYNAMICGRASS_STRING)
					tmp_dclist.push_back(p->get());
			}
			dclist.clear();
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= tmp_dclist.begin(); p != tmp_dclist.end(); p++ )
				dclist.push_back(p->get());
		}
	}


	//创建一个新的点集合，只包含面上的点
	osg::ref_ptr<osg::Vec3Array> geVx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
	{
#if 0			/* 2013-2-26 Galen : 测试发现在前面已经把这个问题解决过了。*/
		/* 
		   如果面上的点跟约束上的点距离很近，面上的这个点就不加进来做三角化，避免在做三角化过程中因为约束上的点跟面上的点距离近而丢失约束上的点并产生空洞的情况出现，
		   by Laker 2013-2-25
		*/
		int nadded=1;
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
		{
			osgUtil::DelaunayConstraint* dc= p->get();
			const osg::Vec3Array* vercon= dynamic_cast<const osg::Vec3Array*>(dc->getVertexArray()); 
			if (vercon)
			{
				for (unsigned int icon=0;icon<vercon->size();icon++)
				{
					osg::Vec3 p1=(*vercon)[icon];
					osg::Vec3 p2=f->pVF;

					double distance = sqrt((p1.x() - p2.x())*(p1.x() - p2.x()) + (p1.y() - p2.y())*(p1.y() - p2.y()));
					if( distance < MINDISTTOMERGEVERTEX)
					{
						nadded = 0;
						break;
					}
				}
			}
			if (!nadded)
				break;
		}

		if (nadded)
#endif
			geVx->push_back(f->pVF);
	}

	osg::ref_ptr<osg::Vec3Array> dcoords = TerrainProcessingUtils::uniqueCoords(geVx.get());		//( coords.get());
	osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(dcoords.get());
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
		dt->addInputConstraint( p->get() );

	dt->triangulate();
	
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		dt->removeInternalTriangles( p->get() );

		bool isCulture = false;
		std::string name = p->get()->getName();
		for (int i=0; i<culture_cc; i++)
		{
			if (name==cultures[i])
			{
				isCulture = true;
				break;
			}
		}
	
		if ( isCulture )
		{
			// Adjust Z coorinates of the terrain matching X,Y in features to the feature Z
			osg::Vec3Array *dcCoords = dynamic_cast<osg::Vec3Array *>((*p)->getVertexArray());
			for( osg::Vec3Array::iterator q = dcCoords->begin(); q != dcCoords->end(); q++ )
			{
				for( osg::Vec3Array::iterator r = dcoords->begin(); r != dcoords->end(); r++ )
				{
					if( (*q)[0] == (*r)[0] && (*q)[1] == (*r)[1] )
						(*r)[2] = (*q)[2];
				}
			}
		}
	} 


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//第四步，1、增加原来去掉的裙边的点，并确定所有点的正确序号

	/* 把三角化后的所有点集合都存在validVertices里面 */
	/* 并把validVertices里面所有点，由投影点映射到空间点 */
	j = 0;
	validVertices.clear();	validElement.clear();
	for( osg::Vec3Array::iterator r = dcoords->begin(); r != dcoords->end(); r++ )
	{
		FindValidVertex vVert;
		osg::Vec3 pSpace = getSpaceVertex(*r);
		vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
		validVertices.push_back(vVert);
	}

#if 0		//如果德罗尼计算正确的话，这段程序可以不要。
	// 如果是去除文化信息内部的无效面，那么要把地形内部与文化信息相交的点都存为无效点。
	for( std::vector<osg::Vec3>::iterator p= innerVerticesInterWithCulture.begin(); p != innerVerticesInterWithCulture.end(); p++ )
	{
		invalidVertices.push_back( *p );
	}
#endif

	//找出所有的有效面。存放在 validElement 里面
	osg::DrawElementsUInt *validdui = dt->getTriangles();
	for(unsigned int i=0; i<validdui->getNumIndices(); i+=3)
	{
		unsigned int a,b,c;
		a = validdui->getElement(i);
		b = validdui->getElement(i+1);
		c = validdui->getElement(i+2);

		if(	(checkTriOnBorder(&validVertices[a], &validVertices[b], &validVertices[c])) &&
				(checkTriOnDclist(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF)) )
		{
			osg::Vec3 normal = getNormal(validVertices[a].pVF, validVertices[b].pVF, validVertices[c].pVF);
			if(normal.z() > 0.0)
			{
				validElement.push_back(a);
				validElement.push_back(b);
				validElement.push_back(c);
			}else{
				validElement.push_back(a);
				validElement.push_back(c);
				validElement.push_back(b);
			}
		}
	}

	//创建一个新的点阵列，加入所有有效面
	osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
		validCoords->push_back(pp->pV);
	
	//把只包含合法三角形的validElement里面的元素加进来
	osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
	for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
		dui->push_back(*p);

	if (dui->empty())	return NULL;

#ifdef ADD_SKIRT
	/* 创建裙边 */
	CreateSkirt(validCoords, dui, 1);
#endif

	/* 创建最终完成切割后的几何体。*/
	osg::ref_ptr<osg::Geometry > newGeometry = new osg::Geometry;   
	newGeometry->setVertexArray( validCoords.get() );
	// Restore texture coordinates
	osg::ref_ptr<osg::Vec2Array> tc = tcp.getTexCoords( *(validCoords.get()) );
	newGeometry->setTexCoordArray( 0, tc.get() );
	// Restore normal 
	osg::ref_ptr<osg::Vec3Array> nc = tcp.getNormals( *(validCoords.get()) );
	newGeometry->setNormalArray( nc.get() );
	newGeometry->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	// Add an overall white color
	osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
	color->push_back(osg::Vec4(1,1,1,1));
	newGeometry->setColorArray(color.get());
	newGeometry->setColorBinding( osg::Geometry::BIND_OVERALL );
	newGeometry->addPrimitiveSet( dui );

#ifdef ADD_SKIRT
	for( std::vector< osg::ref_ptr<osg::DrawElementsUInt> >::iterator p = cultSkirtDrawElements.begin(); p != cultSkirtDrawElements.end(); p++ )
		newGeometry->addPrimitiveSet( *p );
#endif


	/* 如果有细节纹理，就加细节纹理 */
	if(dcDetail.size() > 0)
	{	
		CreateFaceOfDetailTxtWhenHasCulture(newGeometry);
		return NULL;
	}else{
		newGeometry->setStateSet( g->getStateSet() );
		return newGeometry;
	}
}



void processGeode( const GeodeInfo& gi, const GeodeList &features )
{
	if (features.size()<1) return;

	// Set the acceptable fit tolerance, nominally 0.5 units at a distance of
	//   1000 unitss, but scales linearly beyond that.
	double _et = gi._distance / 2000.; // or: distance / 1000. * .5;
	osg::Geode& geode = *(gi._geode.get());
	
	for( unsigned int i = 0; i < geode.getNumDrawables(); i++ )
	{
		osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(geode.getDrawable(i));
		if( geometry.valid() )
		{
			osg::ref_ptr<osg::Geometry> newGeometry = processGeometry( geometry.get(), features );
			if( newGeometry.valid() )
			{
				geode.replaceDrawable( geometry.get(), newGeometry.get() );
			}
			else
			{
				geode.removeDrawable( geometry.get() );
			}
		}
	}

}

/* 增加细节纹理控制层 */
void ChangeCoordSystem(osg::Geode* gd, const osg::Matrix mt);
void AddDetailedTxtIntoTerrain(osg::MatrixTransform* mt, bool bGroup)
{
	/* 把路单独拿出来，*/
	std::vector<osg::ref_ptr<osg::Geode>> road_gdlist;
	std::vector<osg::ref_ptr<osg::Geode>> noRoad_gdlist;
	for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= detailed_gdlist.begin(); p != detailed_gdlist.end(); p++ )
	{
		if(p->get()->getName() == "RoadDay") road_gdlist.push_back(*p);
		else noRoad_gdlist.push_back(*p);
	}
	detailed_gdlist.clear();
	if(noRoad_gdlist.size() > 0)
	{
		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= noRoad_gdlist.begin(); p != noRoad_gdlist.end(); p++ )
			detailed_gdlist.push_back(*p);
	}
	
	/* 如果有路，先处理路。把路也做成多重纹理控制，白天显示原有纹理，夜晚显示发光纹理。*/
	if((mt != NULL)&&(road_gdlist.size()>0))
	{
		osg::ref_ptr<osg::Switch> sw = new osg::Switch;
		sw->setAllChildrenOn();
		sw->setName("IVERoad");

		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= road_gdlist.begin(); p != road_gdlist.end(); p++ )
		{
			if(bGroup == false)	ChangeCoordSystem(p->get(), osg::Matrix::inverse(mt->getMatrix()));
			sw->addChild(p->get());
		}
		road_gdlist.clear();
		mt->addChild(sw.get());
	}


	/* 把湖也单独拿出来，*/
	std::vector<osg::ref_ptr<osg::Geode>> lake_gdlist;
	std::vector<osg::ref_ptr<osg::Geode>> noLake_gdlist;
	for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= detailed_gdlist.begin(); p != detailed_gdlist.end(); p++ )
	{
		if(p->get()->getName().substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING) lake_gdlist.push_back(*p);
		else noLake_gdlist.push_back(*p);
	}
	detailed_gdlist.clear();
	if(noLake_gdlist.size() > 0)
	{
		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= noLake_gdlist.begin(); p != noLake_gdlist.end(); p++ )
			detailed_gdlist.push_back(*p);
	}
	
	/* 如果有湖，先处理湖。把湖。*/
#if 1
	/* IVE文件里面暂时不加入湖面 */
	if(0)	//((mt != NULL)&&(lake_gdlist.size()>0) && (CurrentTerrainLodLevel <= 1))
#else
	if((mt != NULL)&&(lake_gdlist.size()>0))
#endif
	{
		osg::ref_ptr<osgFX::MultiTextureControl> mtc = new osgFX::MultiTextureControl;
		mtc->setTextureWeight(0,0.0f);
		mtc->setTextureWeight(1,1.0f);
		mtc->setName("IVELake");

		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= lake_gdlist.begin(); p != lake_gdlist.end(); p++ )
		{
			if(bGroup == false)	ChangeCoordSystem(p->get(), osg::Matrix::inverse(mt->getMatrix()));
			mtc->addChild(p->get());
		}
		lake_gdlist.clear();
		mt->addChild(mtc.get());
	}

	/* 再处理细节纹理 */
	if((mt != NULL)&&(detailed_gdlist.size()>0))
	{
		osg::ref_ptr<osgFX::MultiTextureControl> mtc = new osgFX::MultiTextureControl;
		mtc->setTextureWeight(0,1.0f);
		mtc->setTextureWeight(1,0.0f);
		mtc->setName("DETAILED_TEXTURE");

		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= detailed_gdlist.begin(); p != detailed_gdlist.end(); p++ )
		{
			if(bGroup == false)	ChangeCoordSystem(p->get(), osg::Matrix::inverse(mt->getMatrix()));
			mtc->addChild(p->get());
		}
		detailed_gdlist.clear();
		mt->addChild(mtc.get());
	}
}


void AddDynamicGrassIntoTerrain(osg::MatrixTransform* mt, bool bGroup)
{
	if((mt != NULL)&&(dynamicGrass_gdlist.size()>0))
	{
		osg::ref_ptr<osg::Group> dyg = new osg::Group;
		dyg->setName("DynamicGrass");

		for( std::vector<osg::ref_ptr<osg::Geode> >::iterator p= dynamicGrass_gdlist.begin(); p != dynamicGrass_gdlist.end(); p++ )
		{
			if(bGroup == false)	ChangeCoordSystem(p->get(), osg::Matrix::inverse(mt->getMatrix()));
			dyg->addChild(p->get());
		}
		dynamicGrass_gdlist.clear();
		mt->addChild(dyg.get());
	}
}


void process( osg::Node* terrain, bool bGroup )
{
	FindFeatures ff;
	terrain->accept( ff );
	
	const GeodeList &features = ff.getFeatures();
	//printf("Features: %d\n", features.size() );
	const GeodeList &nonFeatures = ff.getNonFeatures();
	//printf("Non-Features: %d\n", nonFeatures.size() );
	
	for( GeodeList::const_iterator p = nonFeatures.begin(); p != nonFeatures.end(); p++)
	{
		processGeode( *p, features );
		if (bGroup)
		{
			osg::Node::ParentList parentList =  (*p)._geode.get()->getParents();
			for(osg::Node::ParentList::iterator itr=parentList.begin(); itr!=parentList.end(); ++itr)
			{
				osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(*itr);
				if (mt!=NULL)
				{
					AddDetailedTxtIntoTerrain(mt, true);
					AddDynamicGrassIntoTerrain(mt, true);
					break;
				}
			}
		}
	}
}

void ChangeCoordSystem(osg::Geode* gd, const osg::Matrix mt)
{
	osg::Geode::DrawableList dl = gd->getDrawableList();
	for(osg::Geode::DrawableList::iterator iter = dl.begin(); iter != dl.end(); iter++)
	{
		osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry*> ((*iter).get());
		osg::Vec3Array *vx = dynamic_cast<osg::Vec3Array*>(gm->getVertexArray());
		for(osg::Vec3Array::iterator iter = vx->begin(); iter != vx->end(); iter++)
		{
			osg::Vec3  vo = *iter * mt;
			iter->set(vo.x(), vo.y(), vo.z()); 
		}
	}
}

int UpdateIVEByBTGCSN(osg::ref_ptr<osg::Group> btg, osg::ref_ptr<osg::Node> node)
{
	osg::ref_ptr<osg::Group> root = new osg::Group;
	for(unsigned m=0; m<btg->getNumChildren(); ++m)
	{
		osg::Geode* cgd = dynamic_cast<osg::Geode*>(btg->getChild(m));
		root->addChild(cgd);
	}

	osg::CoordinateSystemNode* csn = findTopMostNodeOfType<osg::CoordinateSystemNode>(node.get());
	if (csn)
	{
		for(unsigned int i=0; i<csn->getNumChildren();++i)
		{
			osg::PagedLOD* plod = dynamic_cast<osg::PagedLOD*>(csn->getChild(i));
			if (plod!=NULL)
			{
				for(unsigned int j=0; j<plod->getNumChildren();++j)
				{
					osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(plod->getChild(j));
					for(unsigned int k=0; k<mt->getNumChildren();++k)
					{
						osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
						
						//直接放入root里面
						root->addChild(gd);

						//所有结点换成世界坐标系 
						ChangeCoordSystem(gd, mt->getMatrix());

						/* 做德罗尼切割*/
						process(root, false);

						/* 所有结点还原为局部坐标系 */
						ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
					}

					AddDetailedTxtIntoTerrain(mt, false);
					AddDynamicGrassIntoTerrain(mt, false);
				}
			}
		}
	}

	return 1;
}


int UpdateIVEByBTGPageLod(osg::ref_ptr<osg::Group> btg, osg::ref_ptr<osg::Node> node)
{
	osg::ref_ptr<osg::Group> root = new osg::Group;
	for(unsigned m=0; m<btg->getNumChildren(); ++m)
	{
		osg::Geode* cgd = dynamic_cast<osg::Geode*>(btg->getChild(m));
		root->addChild(cgd);
	}

	osg::PagedLOD* plod = findTopMostNodeOfType<osg::PagedLOD>(node.get());
	if (plod)
	{
		for(unsigned int j=0; j<plod->getNumChildren();++j)
		{
			osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(plod->getChild(j));
			for(unsigned int k=0; k<mt->getNumChildren();++k)
			{
				osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
				
				//直接放入root里面
				root->addChild(gd);

				//所有结点换成世界坐标系 
				ChangeCoordSystem(gd, mt->getMatrix());

				/* 做德罗尼切割*/
				process(root, false);

				/* 所有结点还原为局部坐标系 */
				ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
			}
			AddDetailedTxtIntoTerrain(mt, false);
			AddDynamicGrassIntoTerrain(mt, false);
		}
	}
	return 1;
}

/* 本地坐标系还原为世界坐标系。（里面可能包含纹理控制层）*/
void RestoreToWorldCoordinate(osg::MatrixTransform* mt)
{
	if(mt)
	{
		for(unsigned int k=0; k<mt->getNumChildren();++k)
		{
			osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
			if(gd != NULL)
			{		
				/* 所有结点还原为局部坐标系 */
				ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
			}else{
				osgFX::MultiTextureControl* mtc = dynamic_cast<osgFX::MultiTextureControl*>(mt->getChild(k));
				if(mtc != NULL)
				{
					for(unsigned int i = 0; i < mtc->getNumChildren(); i++)
					{
						osg::Geode* gd = dynamic_cast<osg::Geode*>(mtc->getChild(i));
						if(gd != NULL)
						{		
							/* 所有结点还原为局部坐标系 */
							ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
						}
					}
				}else{
					osg::Group *gp = dynamic_cast<osg::Group*>(mt->getChild(k));
					if(gp != NULL)
					{
						for(unsigned int i = 0; i < gp->getNumChildren(); i++)
						{
							osg::Geode* gd = dynamic_cast<osg::Geode*>(gp->getChild(i));
							if(gd != NULL)
							{		
								/* 所有结点还原为局部坐标系 */
								ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
							}
						}
					}
				}
			}
		}
	}
}


int UpdateIVEByBTGGroup(osg::ref_ptr<osg::Group> btg, osg::ref_ptr<osg::Node> node)
{
	osg::ref_ptr<osg::Group> root = new osg::Group;
	for(unsigned m=0; m<btg->getNumChildren(); ++m)
	{
		osg::Geode* cgd = dynamic_cast<osg::Geode*>(btg->getChild(m));
		root->addChild(cgd);
	}

	osg::Group* gp = findTopMostNodeOfType<osg::Group>(node.get());
	if (gp)
	{
		for(unsigned int i=0; i<gp->getNumChildren();++i)
		{
			osg::PagedLOD* plod = dynamic_cast<osg::PagedLOD*>(gp->getChild(i));
			if(plod)
			{
				for(unsigned int j=0; j<plod->getNumChildren();++j)
				{
					osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(plod->getChild(j));
					for(unsigned int k=0; k<mt->getNumChildren();++k)
					{
						osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
						
						//所有结点换成世界坐标系 
						ChangeCoordSystem(gd, mt->getMatrix());
	
						//直接放入root里面
						root->addChild(gd);
					}
				}
			}else{
				osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(gp->getChild(i));
				if(mt)
				{
					for(unsigned int k=0; k<mt->getNumChildren();++k)
					{
						osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
						
						//所有结点换成世界坐标系 
						ChangeCoordSystem(gd, mt->getMatrix());

						//直接放入root里面
						root->addChild(gd);
	
					}
				}
			}
		}

		/* 做德罗尼切割*/
		process(root, true);
	
		for(unsigned int i=0; i<gp->getNumChildren();++i)
		{
			osg::PagedLOD* plod = dynamic_cast<osg::PagedLOD*>(gp->getChild(i));
			if(plod)
			{
				for(unsigned int j=0; j<plod->getNumChildren();++j)
				{
					osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(plod->getChild(j));
					RestoreToWorldCoordinate(mt);
				}
			}else{
				osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(gp->getChild(i));
				RestoreToWorldCoordinate(mt);
			}
		}
	}
	return 1;
}

/* 获取当前IVE地景文件的LOD级别 */
void GetCurrentTerrainLodLevel(char *fnIve)
{
	/* 检查输入的地形文件是不是顶级地形文件，并得到不带路径的简单文件名fn. */
	char fn[64];
	int j = 0;
	IsTopTerrain = 1;
	for(unsigned int i=0; i<strlen(fnIve); i++)
	{
		if((fnIve[i] == '\\') || (fnIve[i] == '/'))  { j = 0; }	else { fn[j++] = fnIve[i]; }
	}
	fn[j++] = 0;
	for(unsigned int i=0; i<strlen(fn); i++)
	{
		if(fn[i] == '_') {	IsTopTerrain = 0; break;	}
	}
	
	/* 按下划线将它们分开。*/
	splitMaterialName(fn);
	if(materialParts.size() > 2)
	{
		char levelStr[2]; levelStr[0] = materialParts[1][1]; levelStr[1] = 0;
		CurrentTerrainLodLevel = atoi(levelStr);
	}
	
}

int UpdateIVEByBTG(osg::ref_ptr<osg::Group> btg, char *fnIve)
{
	/* 切割地景前，复位全局通讯变量。*/
	CurrentTerrainIsCut = 0;
	UpdateTheElevationInTerrain = 0;
	allCultures = btg;

	/* 获取当前IVE地景文件的LOD级别 */
	GetCurrentTerrainLodLevel(fnIve);

	/* 初始化纹理信息列表 */
	detailedtxt.clear();
	for(unsigned m=0; m<btg->getNumChildren(); ++m)
	{
		osg::Geode* cgd = dynamic_cast<osg::Geode*>(btg->getChild(m));
		if (cgd!=NULL)
		{
			std::string fname = cgd->getName();
			if(fname.substr(0, strlen(DEFINE_DETAIL_STRING)) == DEFINE_DETAIL_STRING)
			{
				detailedtxt.push_back(fname);
			}
		}
	}
	detailedtxt_cc = detailedtxt.size();

	osg::ref_ptr<osg::Node> te = osgDB::readNodeFile(fnIve);
	osg::CoordinateSystemNode* csn = dynamic_cast<osg::CoordinateSystemNode*> (te.get());
	if(csn)
	{
		UpdateIVEByBTGCSN(btg, csn);
		if (CurrentTerrainIsCut==1)
			osgDB::writeNodeFile(*csn, fnIve);
	}else {
		osg::PagedLOD* plod = dynamic_cast<osg::PagedLOD*> (te.get());
		if(plod)
		{
			UpdateIVEByBTGPageLod(btg, plod);
			if (CurrentTerrainIsCut==1)
				osgDB::writeNodeFile(*plod, fnIve);
		}else{
			osg::Group* gp = dynamic_cast<osg::Group*> (te.get());
			if(gp)
			{
				UpdateIVEByBTGGroup(btg, gp);
				if (CurrentTerrainIsCut==1)
					osgDB::writeNodeFile(*gp, fnIve);
			}
		}
	}
	return 1;
}
