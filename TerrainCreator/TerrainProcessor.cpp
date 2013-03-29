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
/* ��һ���ַ��������»��߻������߷ֳɸ������֡�*/
extern std::vector<std::string> materialParts;
extern int splitMaterialName(std::string matname);
int CurrentTerrainIsCut = 0;
int IsTopTerrain;
int CurrentTerrainLodLevel = 0;
int UpdateTheElevationInTerrain = 0;	//��������˵ؾ�����ĳ����ĸ̡߳�
osg::ref_ptr<osg::Group> allCultures;
double elevationOfLake = -9999.99;		//��¼��ǰ�ĺ��溣�θ߶ȡ�
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

/*  �����ߵ����ƣ����»��߻���������Ϊ�ֶα�־�������Ʒ�Ϊ���Σ�*/
std::vector<std::string> NameParts;
void SplitName(std::string matname)
{
	char matpart[128], mat[256];
	int j = 0;
	
	/* ��matname���»�����Ϊ�ָ������ֳ�����С�Σ�������nameParts���� */
	sprintf(mat, "%s", matname.c_str());
	NameParts.clear();
	for(unsigned int i = 0; i<strlen(mat); i++)
	{
		if((mat[i] != '_')&&(mat[i] != '-')) { matpart[j++] = mat[i];	}
		if((mat[i] == '_')||(mat[i] == '-')) { matpart[j] = 0;	NameParts.push_back(matpart);	j = 0; }
	}
	matpart[j] = 0;	NameParts.push_back(matpart);	j = 0;	
}

/* ��ȡ������ķ��ߡ�*/
osg::Vec3 getNormal(osg::Vec3 v0, osg::Vec3 v1, osg::Vec3 v2)
{
	osg::Vec3 normal = cross(v1 - v0, v2 - v0);
	return normal;
}

/* �õ�һ���߶γ���ֵ */
double distanceV(osg::Vec3 v0, osg::Vec3 v1)
{
	return sqrt( ((v0.x() - v1.x()) * (v0.x() - v1.x())) +  ((v0.y() - v1.y()) * (v0.y() - v1.y())) +  ((v0.z() - v1.z()) * (v0.z() - v1.z())) );
}

/* ���һ���ַ����Ƿ���ϸ�������ַ����� */
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

/* ��һ��Node�����ҵ������ض����� */
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

/* ����Node�������Ļ���Ϣ����ϸ�������node���֡�������_features���棬����ı�����_nonFeatures���� */
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

/* ������Լ�� */
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
bool AllBorderVertexIntersection;				//���б߽�ĵ㶼���Ļ���Ϣ�ཻ��
bool AllInnerVertexIntersection;				//�����ڲ��ĵ㶼���Ļ���Ϣ�ཻ��
bool AllInnerVertexNoIntersection;				//�����ڲ��ĵ㶼���Ļ���Ϣ���ཻ��
int IntersectCt;												//��¼�ؾ�����Ļ���Ϣ�ж��ٸ������ȫ�ֱ�����
int subIntersectCt;											//
bool needDoSkirt;												//����һ����Χ���Ƿ��������ȹ�ߡ�
double skirtLength;											//ȹ�ߵĸ߶� m
osg::BoundingSphere terrainBS;					//�ؾ���İ�Χ��
std::vector< osg::ref_ptr<CorrectionSpaceMesh> > validCSM;

std::vector<ProjectionMap> totalVertices;
std::vector<ProjectionMap> borderPMs;
std::vector<ProjectionMap> featurePMs;
std::vector<ProjectionMap> IntersectVertices;		//��¼���еĽ��㣬�����һ�������ͶӰ�㣬ͬǰ��ĳ������ͶӰ����ͬ����ֱ��ʹ���Ǹ�����Ŀռ����Ϣ���������¼��������㡣
std::vector<ProjectionMap> roadTexCoordSave;		//��¼����·�������������

#ifdef ADD_SKIRT
void CreateSkirt(osg::ref_ptr<osg::Vec3Array> validCoords, osg::ref_ptr<osg::DrawElementsUInt> dui, int mode);
#endif

/* �߶����ݽṹ */
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
/* �������ݽṹ���涥�㵽���߶εľ��롣*/
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
/* �������ݽṹ�����Ѿ��߶λ��ĺ����Ļ���Ϣ */
typedef struct _LAKEEDGEITEM_
{
	int numOfVertex;
	osg::Vec3 firstVertex;
}LakeEdgeItem;
std::vector<LakeEdgeItem> allEdgedLakes;
/* Add by Galen 2013-02-22 End */

FindValidVertex conner[4];
FindValidVertex connerEx[4];						//����conner�ĸ�����ã��������ĸ�����΢����һȦ�ĵ㡣
LineSeg border[4];
std::vector<LineSeg> borderIntersectLS;				//��¼���еı߶Σ�ԭʼ�߶Σ��������Ļ���Ϣ�ཻ(�����˵㶼�ཻ)�����ཻ(ֻ��һ���˵��ཻ)�Ͳ��ཻ(�����˵㶼���ཻ)�����б߶Σ�
std::vector<LineSeg> borderHalfIntersectLS;		//��¼���еİ��ཻ(ֻ��һ���˵��ཻ)�߶Σ�
std::vector<LineSeg> validBorderSegment;			//��¼������Ч�ı߶Σ����������˲����֮�����Ч�ߣ�

int f_idx;

osg::ref_ptr<osg::Vec3Array> FindValidVertexInDC(osg::ref_ptr<ArealConstraint> dc);
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToTranglesFace(osg::ref_ptr<osg::Vec3Array> gvx);

/* ���һ�����Ƿ���FindValidVertex�㼯���� */
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



/* ����conner������connerEx��*/
void calcConnerEx()
{
	Geodesy geo;
	Carton  cat;
	osg::Vec3 v, r;
	double maxLon, maxLat, minLon, minLat;
	maxLon = maxLat = -999999.999;	minLon = minLat = 999999.999;
	ProjectionMap cnEx[4], *pm = &cnEx[0];
	
	/* ������ĸ��ǵľ�γ�ȣ�����ؼ�ֵ��*/
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
	
	/* Ȼ��дConnerEx */
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
	
	/* ���⼸������뵽��Ч�㼯������ */
	invalidVertices.push_back(connerEx[0].pVF);
	invalidVertices.push_back(connerEx[1].pVF);
	invalidVertices.push_back(connerEx[2].pVF);
	invalidVertices.push_back(connerEx[3].pVF);
}

/* ����һ��ͶӰ�㣬�ҵ��õ��Ӧ��3D�ռ�㡣*/
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

/* ����һ��ͶӰ�㣬�ҵ��õ��Ӧ��3D�ռ�㡣*/
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
		/* �߶�ֵ��� */
		geo._elevation +=  totalVertices[i].elev;
		GeodToCart(&geo, &cat);		
		r.x() = cat.x;	r.y() = cat.y;	r.z() = cat.z;
	}
	return r;
}


/* ����һ��Vertex�� ���ת���ɾ�γ������ϵ�󣬽��߶�ֵ��Ϊ0.0�ĵ㣬������õ��ں�ƽ���ӳ��㡣*/
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

/* ����һ��Vertex�� �����ProjectionMap�ṹ*/
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

/* ���ݵ�ǰ���Ƿ��к�����������ǰ�������ĸ߶ȡ�ԭ��͵���csm��߶���ͬ��*/
void TuneElevationOfLakeSurface(osg::ref_ptr<osg::Geode> lakeSurface, osg::ref_ptr<osg::Vec3Array> coords)
{
	if(lakeSurface->getNumDrawables() > 0)
	{
		Geodesy geo;
		Carton  cat;
		osg::Vec3 v;
		
		/* Add by Galen 2013-02-21 ������ߵĵ㡣���ߵ�Ҫ����ƽ������һ�£�����̫ͻأ��*/
		if(vertexNearLake.size() > 0)
		{
			for(std::vector<VertexDistance>::iterator i = vertexNearLake.begin(); i != vertexNearLake.end(); ++i)
			{
				/* ����Ҫ��������㣬��faceVertice�����ҵ�ԭʼ�߶� */
				double currVerHeight;
				for( std::vector<FindValidVertex>::iterator f = faceVertice.begin(); f != faceVertice.end(); f++ )
				{
					if(f->pVF == i->pVF)
					{
						cat.x = f->pV.x(); cat.y = f->pV.y(); cat.z = f->pV.z();
						CartToGeod(&cat, &geo);	currVerHeight = geo._elevation;
						
						/* �������Է�����������ĺ��ߵ�߶ȣ���д�� faceVertice ���� */
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

		/* ����ؾ����ϵĵ� */
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
		
		/* ����ȹ���ϵĵ� */
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
		
		/* ���и߶ȱ仯�ĵ�д�ص�Geometry���档*/
		int index = 0;
		for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
		{
			*r = (terrainVertices[index++].pV);
		}		
	}
}


/* ��֪һ��㣬����ΪFindValidVertex������pV��Ӧ��pVF���������ת��Ϊ��γ��֮�󣬶�Ӧ�߶�Ϊ0���Ǹ��㡣*/
void getFlatVertex(std::vector<FindValidVertex> *plist)
{
	for( std::vector<FindValidVertex>::iterator p= plist->begin(); p != plist->end(); p++ )
	{
		osg::Vec3 pF = getProjectionVertex(p->pV);
		p->pVF = pF; 
	}
}

/* ��֪һ��㣬����ΪFindValidVertex������������ͶӰ��pVF�����һ��osg::Vec3Array */
osg::ref_ptr<osg::Vec3Array> getCoordsFlat(std::vector<FindValidVertex> plist)
{
	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator p= plist.begin(); p != plist.end(); p++ )
	{
		vx->push_back(p->pVF);
	}
	return vx.release();	
}

/* �ѵ�ǰfeature�ĸ���Geometry�������е�Vertex���ĳɾ�γ����������߶�Ϊ0�Ķ�Ӧ�㣬���洢��featurePMs���档 */
void makeGeometryFlat(osg::ref_ptr<osg::Geometry> fg)
{
	osg::ref_ptr<osg::Vec3Array> vx = dynamic_cast<osg::Vec3Array *>(fg->getVertexArray());
	for( osg::Vec3Array::iterator r = vx->begin(); r != vx->end(); r++ )
	{
		/* �����Ļ���Ϣ���ӳ��, ���ں������. */
		ProjectionMap pm;
		getProjectionMapByVertex(&pm, *r);
		featurePMs.push_back(pm);	

		osg::Vec3 pF = pm.pVF;
		*r = pF;
	}
}

/* ����һ��feature�㣬��featurePMs�����ҵ���Ӧ����������Ӧ��浽 totalVertices ���档*/
osg::Vec3 findPMByFeatureVertex(osg::Vec3 vx, int rettype, string pName)
{
	for( std::vector<ProjectionMap>::iterator p = featurePMs.begin(); p != featurePMs.end(); p++ )
	{
		if((vx == p->pV) || (vx == p->pVF))	
		{
			/* ֻ�е��Ļ���Ϣ����Ϊ·��������ϸ���������ʱ�򣬲ŵ����Ļ���Ϣ���ʵ�ʸ߶�ͬ������߶���ͬ��*/
			/* ������������Σ������������Χ�����԰�ԭ�����Ļ���Ϣ��߶ȼ��㡣*/
			if(		(pName.substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING) || 
					(pName.substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) || 
					(isDetailedtxt(pName)) )
			{
				/* ��������һ�δ��룬���¼����Ļ���Ϣ��ĸ߶� 2012-11-25 */
				/* ���㷽����ͬ��ǰ�������������ײ��⣬�޸�featurePMs�����ԭʼֵ��*/
				/* �����ཻ�߶Σ���һ���ھ���lon��γ��lat�Ĵ�ֱ���߶Ρ�������߶���ؾ����ཻ��*/
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
				
				/* ��ȡ�ཻ��� Get intersections */
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
				/* ������ϡ�*/
			}else if (pName.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			/* ���Ļ���Ϣ����Ϊ����ʱ�� Ҳ��ԭ�����Ļ���Ϣ��߶ȼ��㡣 */ 
			{
			}

			totalVertices.push_back(*p);
			if(rettype == 0) 	return p->pVF; else return p->pV;
		}
	}
	return vx;
}

/* ��һ������ ���뵽 IntersectVertices ���� */
void addVertexIntoIntersect(osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, v);
	IntersectVertices.push_back(fvv);	
}

/* �� IntersectVertices ��һ����ͬ�ĵ㣬�Ҳ����ͷ���-1, �ҵ��ͷ����±� */
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



/* ��һ������ ���뵽 totalVertices ���� */
void addVertexIntoTotal(osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, v);
	totalVertices.push_back(fvv);	
}

/* ��һ������ ���뵽 totalVertices ����, �Ȳ������������ totalVertices �����Ƿ������� */
void addVertexIntoTotalFirstSearth(osg::Vec3 v)
{
	Geodesy geo, geo2;
	Carton  cat;
	osg::Vec3 r = v;
	unsigned int minIdx, i, iGet = 0;	//minIdx ������������������Ǹ���֪�߶ȵĶ��㡣
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
/* ��һ������ ���뵽 totalVertices ����, �������Ŀռ�߶�ֵΪ���������� */
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


/* ��һ������Բ��뵽 totalVertices ���� */
void addVertexPairIntoTotal(osg::Vec3 vF, osg::Vec3 v)
{
	ProjectionMap fvv;
	getProjectionMapByVertex(&fvv, vF);
	fvv.pV = v; fvv.pVF = vF;
	totalVertices.push_back(fvv);	
}


/* �����еĶ�����Ϣ���������ε���Ļ���Ϣ�㣬�����е�һ���������棬���Ҽ����ÿ����Ķ�Ӧ��ƽ��㣬��������totalVertices���� */
void getTotalFaceVertex()
{
	/* �ȴ���ؾ����ϵĵ� */
	for( std::vector<FindValidVertex>::iterator p= faceVertice.begin(); p != faceVertice.end(); p++ )
	{
		addVertexIntoTotal(p->pV);
	}
}


/* �����ռ������㣬�������ӳ�䵽��άƽ���ϵ�б�ʺͽض�  */
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

/* ��ֱ������ϵת��Ϊ��γ������ϵ */
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

/* ����γ������ϵת��Ϊֱ������ϵ */
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


/* ����ռ�����ֱ�� AB �� CD ֮����ཻ��� */
// Find intersection between AB-> and CD->
void intersect(  osg::Vec3 A, osg::Vec3 B, 
                 osg::Vec3 C, osg::Vec3 D,
                 osg::Vec3 &I )
{
    double m1, m2;
    double d1, d2;
    bool m1i, m2i;  // m1 is not infinite, m2 is not infinite;

#if 1
	/* �Ȱ�ֱ������ϵת��Ϊ��γ������ϵ */
	osg::Vec3d gA, gB, gC, gD, gI;
	gA = cTog(A);	gB = cTog(B);	gC = cTog(C);	gD = cTog(D);

	/* �������γ�ȸ�ʽ�Ľ��� */
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
		�ռ�ֱ�ߵ�����ʽ�� 
		��������ƽ������ϵ�е�����ʽ�� 
		������ΪA(x1,y1,z1),B��x2,y2,z2)
		��ֱ��AB����Ϊ(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
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
/* ��ȡ�ռ�������A��B֮����ֱ��AB�ϵĵ�����(10)�������ȵĵ㣬������vec_list����  */
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
		/* fabs(delx)��fabs(dely)��ȣ�����������0.0 */
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

/* �ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
void FindOutOddPoints(std::vector<SimEdge> *edges)
{
	/* ��һ�γ���Ϊ��ɾ���߻�������������Ȧ��ʹ���߻����һ�������߻���*/
	std::vector<VerIdx> AllPoints;
	std::vector<VerIdx> OddPoints;
	std::vector<VerPair> OddPointPairs;

	/* ���ȣ�������ظ��ıߣ���ȥ����*/
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
		
	/* ��һ�����ҳ��߼��������еĶ��� */
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
		
		/* �ڶ���Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		OddPoints.clear();
		for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
		{
			if(r->times == 1)
			{
				VerIdx vi; vi.Ver = r->Ver; vi.idx = r->idx; vi.times = r->times;
				OddPoints.push_back(vi);
			}
		}
		
		/* ����һ�����⴦��״�����������ֻ������OddPoints��˵������߻��ǷǷ�յ��߻���*/	
		/* ���������������������ʹ���߻����һ������߻���2012-10-27 */
		if(OddPoints.size() == 2)
		{
			/* ��������������ӳ�һ���Σ�*/
			SimEdge se;	 
			se.A = OddPoints[0].Ver; se.B = OddPoints[1].Ver; se.x = 0; se.y = 0; se.z = 0;
			edges->push_back(se);
			OddPoints.clear();
		}

		/* �������֮����û��ֱ�����ӵ��߶Σ��оͰ�����߶�ȥ���� */
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
			
			/* ��������totalEdges */
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


/* �ҳ���ͷ�㲢��������ͷ�� */
void RemoveOddPoints(std::vector<SimEdge> *edges)
{
	/* ��һ�γ���Ϊ��ɾ���߻�������������Ȧ��ʹ���߻����һ�������߻���*/
	std::vector<VerIdx> AllPoints;
	std::vector<VerIdx> OddPoints;
	std::vector<VerPair> OddPointPairs;
		
	/* ��һ�����ҳ��߼��������еĶ��� */
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
		
		/* �ڶ���Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		OddPoints.clear();
		for( std::vector<VerIdx>::iterator r = AllPoints.begin(); r != AllPoints.end(); r++ )
		{
			if(r->times == 1)
			{
				VerIdx vi; vi.Ver = r->Ver; vi.idx = r->idx; vi.times = r->times;
				OddPoints.push_back(vi);
			}
		}
		
		/* �������֮����û��ֱ�����ӵ��߶Σ��оͰ�����߶�ȥ���� */
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
			
			/* ��������totalEdges */
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



/* �� validBorderSegment ���һ���߻� */
void CreateLineLoopFromValidBorderSegment(osg::ref_ptr<osg::Vec3Array> coords)
{
	std::vector<SimEdge>  totalEdges;	totalEdges.clear();
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		SimEdge se;	 
		se.A = q->pS; se.B = q->pE; se.x = 0; se.y = 0; se.z = 0; se.keep = true;
		if(se.A != se.B)	totalEdges.push_back(se);
	}

	/* �ҵ�����������*/
	FindOutOddPoints(&totalEdges);
	
	/* �п��� validBorderSegment ���ϲ���һ����յĶΣ���ô��ʱ�� totalEdges ����Ŀ�ͻ���0 */
	if(totalEdges.size() == 0)	{		return;	}

	while(totalEdges.size() > 2)
	{
		/* �����е�����һ����LINE_LOOP*/
		for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ ) p->keep = true;
		SimEdge start = totalEdges.front();

CPA_LOOP_0:
		SimEdge e = start;
		osg::Vec3 A = e.A;
		osg::Vec3 B = e.B;
		osg::Vec3 firstA = A;
		bool bEnd = true;	//�����ǰ�����бߴղ���һ���߻���bEnd�ͻ����true�����������˳�do..whileѭ��������ܴճ�һ���߻���bEndʼ��Ϊfalse.
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
	
		/* �������ѭ������bEndΪ�����˳�ѭ����Ļ���˵���õ�������߻������Ǻϸ���߻���2012-10-11 */
		if(e == start)	coords->push_back(firstA);
		else
		{
			/* �� e ����ߣ���ԭ��totalEdges �����õ���Ȼ�����¼��㡣 */
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


			/* ������LINELOOP */
			coords->clear();
			goto CPA_LOOP_0;
		}
		
		/* �ҵ�һ���߻�֮���Ȱ�����߻����������浽 multiDc*/
		if(coords->size() > 3)	
		{
			osg::ref_ptr<osg::Vec3Array> nCords = dynamic_cast<osg::Vec3Array *> (coords->clone(osg::CopyOp::DEEP_COPY_ALL));
			multiDc.push_back(nCords);
		}
		coords->clear();

		/* �� totalEdges ���� keep == false�ı�����һ����Ȼ�����¼��㡣 */
		{
			/* ������һ��false�� */
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
			/* ʣ�µı���ȥ�� ��� */
			RemoveOddPoints(&remain);
			/* ʣ�µı��ٽ�����totalEdges ���¼�����һ���߻�*/
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
	int index;		//���Ļ���Ϣ�ཻ�ı߽�ε���ţ�
	int num;		//�Ļ���Ϣ����α߽��ཻ�˼��Σ�
}SegmentIntersection;

typedef struct _DISTANCEOFSEGMENT_
{
	osg::Vec3	vertex;		//����
	double		distance;	//������ζ˵�ľ��롣
}DistanceOfSegment;

/* ����Ҫ����һ�����������Ļ���Ϣ������һ�α߽��ж�ε��ཻ������ô�죿*/
std::vector<int> removeValidSegmentWhenWithMultiIntersections;
std::vector<FindValidVertex> intsResult;  		//��¼���Ļ���Ϣ��Χ������α߽��ཻ�ĵ㼯�ϡ�
void CheckSegmentsIntersection()
{
	unsigned int i,j,k,m;
	std::vector<SegmentIntersection> allSis;	allSis.clear();
	
	/* ��Ҫ���һ�� removeValidSegmentWhenWithMultiIntersections �����Ƿ�����ж����߶��ཻ��*/
	/* ��Ϊ�ж�removeValidSegmentWhenWithMultiIntersections�������ǣ�*/
	/* �߽磺n-1���ཻ��n�㲻�ཻ��˵�� n-1 �� n ����֮��������߶�ͨ����������ǿ��ܡ�*/
	/* ��ʱ��� n �����Ļ���Ϣ�����غϣ���ô�ͻ������һ��n-1��n���߶Σ�����������ܵ�ԭ���ཻ�߽��߶β����غ� */
	/* ���ǲ��ܱ�֤ n-1 �� n ����֮��ȷʵ���߶�ͨ�������Ի�Ҫ���һ�¡�*/
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

	/* �������ཻ���� validBorderSegment ����ͳ������ */
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
	
	/* �ж�ͳ�ƽ��allSis��size()��validBorderSegment��size()�Ƿ���ȣ�
	   �����ȣ�˵��ÿһ���ཻ��ֻ��һ�����㣬�����������������
	   �������ȣ�˵��ĳЩ�ཻ�����ж�����㡣
	 */
	if( allSis.size() < validBorderSegment.size() )
	{
		for(i=0; i<allSis.size(); i++)
		{
			if(allSis[i].num > 1)
			{
				/* �ҳ�����ཻ�ε������ཻ�㣬��������ý���Ͷζ˵�ľ��룬�ռ�������allDOS���� */
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
							/* �ҳ���һ�߽�ε�ԭʼ�������˵㣬������S_fst��S_end���� */
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
				
				/* IntFlag == 2 ��ʾԭʼ�߽�������˵㶼���Ļ���Ϣ�ཻ */
				/* ��Ҫ�ж�����߽�������˵��Ƿ��ཻ״����Ϊ1�������Ϊ1����¼������Σ��Ա���߰���ȥ����*/		
				if(IntFlag == 2)
				{
					DistanceOfSegment dos;	dos.vertex = S_end; dos.distance = LenOf2V(S_fst, S_end);	allDOS.push_back(dos);
					removeValidSegmentWhenWithMultiIntersections.push_back(allSis[i].index);
				}

				/* ������������������Ԫ�أ���distance��С����������� */
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
				
				/* ������֮�󣬽���Ч�μ��뵽 validBorderSegment ���� 
				   ����������һ�����ȰѼ����������Ч�μ��뵽��ʱ����currBorderSegment����
				   �ڶ�������ԭ�� validBorderSegment �������Ч��ȫ��ɾ����
				   ���������� validBorderSegment �� currBorderSegment ��Ч�������һ�� */

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
				
				/* ���� */
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

	/* ���һ����Ч�߶ε������˵���ͬ��ȥ������߶Ρ�*/
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

/* ������������ */
typedef struct _BorderVertexStatus_
{
	int idx;
	int isValid;	
}BorderVertexStatus;
std::vector<int> mergeVertex;			//����߽�ζ˵�ͽ���ĺϲ�״��

void CheckSpecialStatus()
{
	/* Add by Galen 2013-02-05 */
	/* ��� validBorderSegment ����������ͬ��˵���������Ļ���Ϣ������αߵĽ����غϡ����������Ȼ�����˵��α��ϣ����ǲ�δ��������α߽���ʵ���Ե��и*/
	/* Ӧ������������ͬ��ȥ�������Ұ��Ǹ�����Ҳ�ӽ��㼯�� allBorderVertices ����ȥ����*/
	int hasSameSegment = 0;
	std::vector<osg::Vec3> interpVFs;
	for( std::vector<LineSeg>::iterator q = validBorderSegment.begin(); q != validBorderSegment.end(); ++q )
	{
		/* �ȸ���һ�� */
		std::vector<LineSeg> tmpSegment;
		for( std::vector<LineSeg>::iterator p = validBorderSegment.begin(); p != validBorderSegment.end(); ++p )
			if(p != q)	tmpSegment.push_back(*p);
		
		/* �ٱȽ�, ע�⣺����������Ӧ���ٶ��һ���ж����������ǣ��������ڱ߶������˵����Ļ���Ϣ����ཻ״̬Ӧ������ͬ�ġ��Ժ��ټ� 2013-02-06 */
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
		/* �Ȱ���ͬ�ߴ�validBorderSegment����ȥ�� */
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

		/* �ٴӽ��㼯������ȥ����صĽ��㡣*/
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

	/* ����һ��validBorderSegment����� validE ֵ����������״���£�*/
	/* validBorderSegment[].validE = 0(�ȼ��������)������� pE��Ҳ�ǽ���(����������)����ôҪ����� validBorderSegment[].validE ��Ϊ 1 */
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		if(q->validE == 0)
		{
			int curIdx = q->idxE;
			if(terrainBorderVertices[curIdx].newIdx == 1)
				q->validE = 1;
		}
	}

	/* ������������ֻ���������㣬���ڲ�ͬ�Ķ��ڣ�������������ͬ�Ķη���������ͬ�ıߣ���ô����������һ������һ���ǡ�
	   �����취����ԭ���������ӵ������������� */
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


	/* ������������ֻ���������㣬������һ�����ڣ���ô������ε���β�����˵㣬�ཻ״̬Ӧ����һ�µģ��������һ�µľʹ��ˡ�
	   �����취�����validE״̬��1����ô��validE״̬�ĳ�0. */
	if(	(validBorderSegment.size() == 2) && 
		(validBorderSegment[0].border == validBorderSegment[1].border) &&
		(validBorderSegment[0].validE == 1) &&
		(validBorderSegment[1].validE == 1) )
	{
		validBorderSegment[0].validE = validBorderSegment[1].validE = 0;
	}

	/* ������������������Ч�ཻ�θ���Ϊ0�����ǰ��ཻ�εĸ���ȴ��Ϊ0��˵���ཻ���̫С�ˣ��ɺ��Բ��ơ�
	   �����취�������߽�ĸ����˵���ཻ״̬��
	   ����������ཻ���ڲ�ͬ�ıߣ���������������Ͳ����ڡ�������Ч�߶�Ϊ0�������ཻ�����Ϊ0�������ж��ͱ߽���غϵ������
	   ֻ�������ཻ����ͬһ���ϣ�������������ų�����*/
	if((borderHalfIntersectLS.size() > 0) && ( validBorderSegment.size() == 0))
	{
		/* ���жϽ��������intsResult[]��������������������������������ͬ��ͬ�ı��ཻ����ô������������Ͳ�������*/
		if(!((intsResult.size() >= 2) && (intsResult[0].isValid != intsResult[1].isValid)))
		{
			/* ��ǰ�߽�����ж��ٸ��ཻ��*/
			int intBordVers = 0;
			for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
			{
				if(q->newIdx != 0) intBordVers++;
			}
	
			/* �����ǰ�߽�����ཻ��������� 10 ����*/
			if(intBordVers < 10)
			{
				/* �ѱ߽��ȫ����ཻ��*/
				for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
				{
					q->newIdx = 0;
				}
			}
	
			/* ��������е���Ч�ཻ�߶� */	
			validBorderSegment.clear();
			borderHalfIntersectLS.clear();
			return;
		}
	}

	/* �������������ཻ�εĸ������ڽ���ĸ�����˵���н�����㲻�������Ļ���Ϣһ���ߺͱ߽�����غϡ�
	   �����취�������߽�ĸ����˵���ཻ״̬�� */
	if(borderHalfIntersectLS.size() >  (validBorderSegment.size() + mergeVertex.size()))
	{
		/* ���ʣ����������ཻ�����ڣ��ҷ�����ͬ�ıߣ�˵�����������ཻ�����е���һ���ǵ㡣*/
		/* ����ǵ�������ཻ��û�м�������㣬˵������ǵ������������Ļ���Ϣ���ϡ���ô������ǵ㻹ԭ���ɡ�*/
		std::vector<LineSeg> validHalfIntersectLS;
		for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
		{
			if(r->border != -2)
			{
				validHalfIntersectLS.push_back(*r);
			}
		}
		/* ֻʣ��������Ч�İ��ཻ�� */
		if(validHalfIntersectLS.size() == 2)
		{
			/* ��������Ч���ཻ��Ҫ���� */
			if ((validHalfIntersectLS[0].idxS == validHalfIntersectLS[1].idxE) ||
					(validHalfIntersectLS[0].idxE == validHalfIntersectLS[1].idxS) )
			{
				/* ��������Ч���ཻ�η�����ͬ�ı� */
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
		
		/* �ȿ����Ƿ������ڲ��㶼���Ļ���Ϣ���ཻ������ཻ����˵�����غ�����*/
		if(AllInnerVertexIntersection)
		{
			/* ��ǰ�߽�����ж��ٸ��ཻ��*/
			int intBordVers = 0;
			for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
			{
				if(q->newIdx != 1) intBordVers++;
			}

			/* �����ǰ�߽���в��ཻ��������� 3 ����*/
			if(intBordVers < 3)
			{
				/* �ѱ߽��ȫ����ཻ��*/
				for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
					q->newIdx = 1;
				/* ��������е���Ч�ཻ�߶� */	
				validBorderSegment.clear();
			}
			return;
		}
		
		/* ������һ�����ڲ������Ļ���Ϣ���ཻ����� */
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
			
			/* ȥ����Ч��validBorderSegment�� */
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

	/* ������������һ���߽�ε������˵㶼���Ļ���Ϣ�ཻ���߲��ཻ����������߽����ȴ�ж�ֻ��һ�����㣬���������Ȼ�Ǵ���ġ�����ȴ������!!!.
	   �����취���ҵ�������ཻ�Σ��޸ı��κ�����ཻ��֮����ཻ��Ϣ�� */

	/* ����Ҫ�ҵ������˵��ཻ״̬��ͬ�����Ƕ���ȴ�н�������ж� */
	std::vector<int> bothSegment;
	std::vector<int> diffSegment;
	for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
	{
		if(q->validE != 1) bothSegment.push_back(q->border);
		if(q->validE == 1) diffSegment.push_back(q->border);
	}
	
	/* ����߽�������˵��ǲ�ͬ�ཻ״̬��������� */
	if(diffSegment.size())
	{
		/* ���Ҫ�ҵ������˵��ཻ״̬��ͬ�����Ƕ���ȴ��ż������������д���� */
		std::vector<int> errorSegment;
		for(unsigned int i=0; i<diffSegment.size(); i++)
		{
			int iNum = 0;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				if(q->border == diffSegment[i]) iNum++;
			if(iNum%2 == 0) errorSegment.push_back(diffSegment[i]);
		}
		
		/* Ȼ�����Щ����ν��д��� */
		if(errorSegment.size())
		{
			/* ����errorSegment�ж��ཻ״̬ */
			for(unsigned int i=0; i<errorSegment.size(); i++)
			{
				for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				{
					if(q->border == errorSegment[i])
					{
						/* ����validBorderSegment �ཻ״̬ */
						/* ����terrainBorderVertices��Ӧֵ�� �ཻ״̬ */
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
	
	/* ����߽�������˵�����ͬ�ཻ״̬��������� */
	if(bothSegment.size())
	{
		/* ���Ҫ�ҵ������˵��ཻ״̬��ͬ�����Ƕ���ȴ����������������д����, ����һ���޶��������������ֵ����1��  */
		std::vector<int> errorSegment;
		for(unsigned int i=0; i<bothSegment.size(); i++)
		{
			int iNum = 0;
			for( std::vector<LineSeg>::iterator q= validBorderSegment.begin(); q != validBorderSegment.end(); q++ )
				if(q->border == bothSegment[i]) iNum++;
			if((iNum%2 == 1)&&(iNum > 1)) errorSegment.push_back(bothSegment[i]);
		}

		/* Ȼ�����Щ����ν��д��� */
		if(errorSegment.size())
		{
			/* �ж����а��ཻ�߽�Σ��Ƿ���û�������İ��ཻ�߽�Ρ�*/

			/* �ҳ����а��ཻ������û�н���Ķε������˵㼰״ֵ̬��*/
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

			/* ����errorSegment�ж��ཻ״̬ */
			for(unsigned int i=0; i<errorSegment.size(); i++)
			{
				int currFlag = terrainBorderVertices[errorSegment[i]].newIdx;
				int currIdx  = terrainBorderVertices[errorSegment[i]].idx;
				int minDist  = 999, minIdx = -1;
				
				/* �ҵ��뵱ǰ�����ཻ����̾���İ��ཻ�߽�εĶ˵� */
				for(unsigned int j=0; j<allBVS.size(); j++)
				{
					if(currFlag == allBVS[j].isValid)
					{
						int dist = abs(currIdx - allBVS[j].idx);
						if(dist < minDist) {	minDist = dist; minIdx = allBVS[j].idx ;	}
					}
				}
				
				/* ת��terrainBorderVertices���沿�ֵ���ཻ״̬�� */
				unsigned int m, n;
				if(minIdx <= currIdx) {	m = minIdx;	n = currIdx;	}
				if(minIdx > currIdx) {	m = currIdx + 1;	n = minIdx;	}
				for(unsigned int j=m; j<=n; j++)
				{
					if(currFlag==0)	terrainBorderVertices[j].newIdx = 1;
					if(currFlag==1)	terrainBorderVertices[j].newIdx = 0;
				}
				
				/* ���ı�validBorderSegment������Ӧ�ε�״̬ */
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

/* �����һ�����������˽ǵ��� */
void CheckSpecialIntVertex(osg::Vec3 I)
{
	int iFlag = 0;
	for(unsigned int i=0; i<4; i++)
	{
		double dist = LenOf2V(I, conner[i].pVF);
		if(dist < MINDISTTOMERGEVERTEX)	{ iFlag = 1; break;}
	}
	if(!iFlag) return;
	
	/* �ҳ�����ǵ����ڵİ��ཻ�Ρ�*/
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

	/* ������������������ཻ�ζ���������ǵ㣬���˳��� */
	if(validHalfIntersectLS.size() != 2) return;
#if 0
	/* ��δ������ͼ��֪�ιʡ������ˣ�����ΪʲôҪ�ѡ������ǵ���ཻ�Ρ���ɾ�����ȹص���2012-10-19*/
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


/* �� checkBorderVertices �������һ�����㣬���㣺A:onCulture = 1; B:������������������С��1.0 .*/
/* ����ҵ������ĵ㣬������onCulture = 0 */
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

/* ���кϲ����������ֺ�Ҫ�� faceVertice �����������Ӧ�ĵ�Ҫ�޸�һ�¡����򣬺��������������ǻ���ʱ��ͻ�ò�����ȷ�����*/
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

/* �ж�����һ������������Ч�ཻ����Ϊ0���ཻ��ֻ��һ��������������ǽǵ㣬��ô�����ж�����Ļ���Ϣ�����ؾ����ཻ��*/
bool DetectIntersectStatusOne()
{
	/* �����Ч�ཻ�θ���������0��ֱ�ӷ���false��*/
	if(validBorderSegment.size() > 0)	return false;
	
	/* ͳ���ཻ�����������1��������false */
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
	
	/* �ж����Ψһ���ཻ���ǲ��ǽǵ㡣 �������,����false */
	int iFlag = 0;
	osg::Vec3 I = fvv.pVF;
	for(unsigned int i=0; i<4; i++)
	{
		double dist = LenOf2V(I, conner[i].pVF);
		if(dist < MINDISTTOMERGEVERTEX)	{ iFlag = 1; break;}
	}
	if(!iFlag) return false;
	
	/* �������������㣬���� true */
	return true;
}

/* �����������ɵ��߶��м䣬ͨ���𲽵���ֵ���뷨���ҵ�һ�������һ����С��10�׵Ĳ���㡣*/
osg::Vec3 getIntersectVertexLessTenMeter(osg::Vec3 src, osg::Vec3 tgt)
{
	osg::Vec3 mid;	mid.x() = (src.x() + tgt.x()) / 2.0;	mid.y() = (src.y() + tgt.y()) / 2.0;	mid.z() = (src.z() + tgt.z()) / 2.0;
	double length = LenOf2V(src, mid);
	if(length < 10.0)
		return mid;
	else
		return getIntersectVertexLessTenMeter(src, mid);
}

/* ����߻�pm �� ����gd ���ཻ������������ཻ�㵽�߻����� */
/* ���ص��߻�һ��Ҫ���ں�ƽ��ͶӰ���߻�. */
void addIntersectVertex(osg::ref_ptr<osg::Geode> gTerrain, osg::ref_ptr<osg::Geode> gCulture, osg::ref_ptr<osg::Vec3Array> pm, osg::ref_ptr<osg::Vec3Array> newpm, string pmName)
{
	if (pm->size()<3) return;

	/* ��pm����������ȱ��浽 totalVertices ���� */
	for( osg::Vec3Array::iterator r = pm->begin(); r != pm->end(); r++ )
		*r = findPMByFeatureVertex(*r, 1, pmName);

	std::vector<FindValidVertex> fPerimeter;		//��¼һ��ͶӰ����ƽ�����Χ��
	std::vector<FindValidVertex> fPerimeterNoFlat;	//��¼һ��ԭ������Χ��
	std::vector<LineSeg> intersectLS;				//��¼���ཻ�ı�
	std::vector<LineSeg> intersectLS_NoFlat;		//��¼���ཻ�ı�

	intsResult.clear();
	fPerimeter.clear();
	fPerimeterNoFlat.clear();
	removeValidSegmentWhenWithMultiIntersections.clear();
	multiDc.clear();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �鿴�м����ǵ���Լ����
	int ic;
	int icc=0;
	for(ic=0;ic<4;ic++)
	{
		if(conner[ic].isValid > 0) {conner[ic].isValid = 1;	icc++; }
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* �ҳ�0,3������ϵ��������꣬����arctan������������ˮƽ�ߵļнǣ��������ת���� */
	double alpha = 0;
	if (conner[3].pVF[1]-conner[0].pVF[1] != 0)
		alpha = atan((conner[3].pVF[0]-conner[0].pVF[0]) / (conner[3].pVF[1]-conner[0].pVF[1]));
	osg::Matrix _mr = osg::Matrix::rotate(alpha, osg::Z_AXIS);
	osg::Vec3  r_conner[4];
	for (ic=0;ic<4;ic++)
		r_conner[ic] = conner[ic].pVF * _mr;

	int numFi = 0;

	/* ʹpm�߻���֤��β������ */
	if (!pm->empty() && pm->front()!=pm->back())
		pm->push_back(pm->front());

	/* �������pm�߻����鿴�߻��ϸ�������ؾ���gd���ཻ��ϵ��������ཻ�����ʾ�õ���ӳ�䵽gd���ϡ�  */
	for( osg::Vec3Array::iterator p = pm->begin(); p != pm->end(); p++ )
	{
		/* ����ཻ����ʾ������ܴ�ֱӳ�䵽 gTerrain ���ϡ�	*/
		FindValidVertex fvv;

		/* �ҵ���Ӧ��ͶӰ�� */
		osg::Vec3 pFlat = getProjectionVertex(*p);

		osg::Vec3 isect;
		if( TerrainProcessingUtils::getIntersect( gTerrain, pFlat, isect ) != false )
		{
			fvv.isValid = 1;
		}
		else
		{
			//////////////////////////////////////////////////////////////////////////////////////
			// ��������֮��ľ���̫�󲢿�����Σ�������ȷ�ж��ཻ�������

			if (numFi>0 && fPerimeter[numFi-1].isValid==0)
			{
				osg::Vec3 A, B;

				/* �ҵ���ǰ��������ڵ�֮ǰ��һ����  */
				A = pFlat;
				B = fPerimeter[numFi-1].pV;

				/* ͨ��������ת�����ռ�������ת��Ϊƽ��ֱ������ϵ�ϵ������㣬�Ա�������������ƽ��ֱ������ϵ������������  */
				osg::Vec3 rA   = A * _mr;
				osg::Vec3 rB   = B * _mr;

				/* ����������Ƿ���ͬһ������, �����ͬһ���ޣ��򲻴��� */
				bool bProcess = false;
				for (ic=0;ic<4;ic++)
				{
					/* ԭ���ƶ��� rorg �㣬�����ĸ��ǵ�ֱ���Ϊԭ�� */
					osg::Vec3 roA  = rA - r_conner[ic];
					osg::Vec3 roB  = rB - r_conner[ic];

					/* ���ֱ������ϵ���������x������Ų�ͬ ����y������Ų�ͬ������������϶�����ͬһ��������  */
					if ( (roA[0]>0 && roB[0]<0) || (roA[0]<0 && roB[0]>0) || (roA[1]>0 && roB[1]<0) || (roA[1]<0 && roB[1]>0) )
						bProcess = true;

					if (bProcess) break;
				}

				/* ������㲻��ͬһ������, ����A,B��֮����в�ֵ */
				if (bProcess)
				{
					/* ��ֵ��ƽ������10����  */
					std::vector<osg::Vec3> vec_ab;
					interpolation(A, B, vec_ab, INSERT_SEGMENTS);

					/* ����10������ؾ��� gTerrain ���ཻ�����nB��ʾ��ʼ�ཻ�����ţ�nE��ʾ�����ཻ������  */
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

					/* �����ֵ����vec_ab��ؾ���gd���ཻ��  */
					if (nB>=0 || nE>=0)
					{
						/* ��������ཻ�����е��е�  */
						osg::Vec3 iVec;
						if (nB>=0 && nE>=0)
							iVec = (vec_ab[nB]+vec_ab[nE])/2;
						else if (nB>=0)
							iVec = vec_ab[nB];

						/* ������е���뵽��Χ��fPerimeter��  */
						FindValidVertex fvv_n, fvv_nF;
						fvv_nF.pV = iVec;  fvv_nF.idx = numFi; fvv_nF.isValid = 1;
						fPerimeter.push_back(fvv_nF);

						/* �����ԭ��ֵ������Ĳ�ֵ�� */
						A = *p;
						B = fPerimeterNoFlat[numFi-1].pV;
						/* ��ֵ��ƽ������10����  */
						std::vector<osg::Vec3> vec_ab;
						interpolation(A, B, vec_ab, INSERT_SEGMENTS);
						if (nB>=0 && nE>=0)
						{
							if(nE > (int)(vec_ab.size() - 1)) nE = vec_ab.size() - 1;
							iVec = (vec_ab[nB]+vec_ab[nE])/2;
						}
						else if (nB>=0)
							iVec = vec_ab[nB];
						/* ������е���뵽��Χ��fPerimeterNoFlat��  */
						fvv_n.pV = iVec;  fvv_n.idx = numFi; fvv_n.isValid = 1; numFi++;
						fPerimeterNoFlat.push_back(fvv_n);

						/* �Ѷ���Լ��뵽totalVertices���� */
						addVertexPairIntoTotal(fvv_nF.pV, fvv_n.pV);
					}
				}
			}
			//////////////////////////////////////////////////////////////////////////////////////

			fvv.isValid = 0;
		}
		
		/* ����ǰ����뵽��Χ��fPerimeter�� */
		fvv.pV = pFlat;  fvv.idx = numFi;
		fPerimeter.push_back(fvv);

		fvv.pV = *p;  fvv.idx = numFi; numFi++;
		fPerimeterNoFlat.push_back(fvv);
	}
	
	/* ������֮���ҳ��ཻ�߶Σ�isValid��0��1�ʹ�1��0�任�ĵط������ཻ�߶� */
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

	/* �ҳ�ƽ����Ӧ��ά�ռ����ཻ�߶� */
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


	/* ����н���Լ����������û���ཻ�߶Σ����ҵ��ο���ڲ��㶼�����Ļ���Ϣ�ཻ����ôֻ��һ�����������Ļ���Ϣ����ο�ֻ��һ�����㣬�����Ǹ��ǵ㡣*/
	if((icc > 0) && (intersectLS.size() == 0) && (AllInnerVertexNoIntersection == true)) 
	{
		/* �ҵ��Ǹ��ཻ�Ľǵ㣬������onCulture״̬�Ļ�����ʹ����㲻�����ൽ invalidVertices �������档*/
		for( std::vector<FindBorderVertex>::iterator p= checkBorderVertices.begin(); p != checkBorderVertices.end(); p++ )
		{
			if(p->onCulture == 1)	p->onCulture = 0;
		}
		
		/* �˳� */
		pm->clear();
		newpm->clear();
		return;
	}

	/* ����ĸ��Ƕ���Լ����������û���ཻ�߶Σ����ҵ��ο�����е��ڲ��㶼���Ļ���Ϣ�ཻ����ô����˵���ؾ�ȫ�����Ļ���Ϣ�ڲ�*/
	/* ��������ؿ���ĸ������һ��LINE_LOOP��� */
	/* �������Ļ���Ϣ�Ǻ��棨���棩�����ҵؾ���������ڲ��㶼������Ļ���Ϣ�ཻ����ô�ͱ�ʾ �����ؾ��涼�Ǻ��棨���棩----- �������� */
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
		
		/* ���� checkBorderVertices �����е��ֵཻ����ʾ�������еĵ㶼���Ļ���Ϣ�ཻ��*/
		if(	(pmName.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING) || 
				(pmName.substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) )
		{
			/* ����ȷ������߻� */
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

	/* ������еı߽�ζ����Ļ���Ϣû�н��㣬�Ǿ�ֱ���˳��������ٿ�����һ�����������Ļ���Ϣ����ȫ�ڵؾ����ڡ�*/
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

	/* �����Щ�ཻ�߶�����ߵ��ཻ�������������� intsResult ���� */
	validBorderSegment.clear();
	mergeVertex.clear();
	for(j=0; j<(int)(intersectLS.size()); j++)
	{
		// ��dif��֤һ���߶θ�������һ����һ������		
		bool bInter = false;
		int dif = 0;
		while (bInter == false)
		{
			for(i=0; i<4; i++)			//i = intersectLS[j].border;		
			{
				osg::Vec3 I;

#if 1
				/* Ҫʹ I ����intersectLS���ϣ���Ϊ�����Ļ���Ϣ���ߡ�*/
				intersect(intersectLS[j].pS, intersectLS[j].pE, border[i].pS, border[i].pE, I);
#else
				/* Ҫʹ I ����border���ϣ���Ϊ���ǵ��α߽���ߡ�*/
				intersect(border[i].pS, border[i].pE, intersectLS[j].pS, intersectLS[j].pE, I);
#endif

				//�ж�I�Ƿ����߶�intersectLS[j]��
				if(   (((I[0] >= intersectLS[j].pS[0]-dif)&&(I[0] <= intersectLS[j].pE[0]+dif)) || ((I[0] <= intersectLS[j].pS[0]+dif)&&(I[0] >= intersectLS[j].pE[0]-dif)))
					&&(((I[1] >= intersectLS[j].pS[1]-dif)&&(I[1] <= intersectLS[j].pE[1]+dif)) || ((I[1] <= intersectLS[j].pS[1]+dif)&&(I[1] >= intersectLS[j].pE[1]-dif)))
					&&(((I[0] >= border[i].pS[0]-dif)&&(I[0] <= border[i].pE[0]+dif)) || ((I[0] <= border[i].pS[0]+dif)&&(I[0] >= border[i].pE[0]-dif)))
					&&(((I[1] >= border[i].pS[1]-dif)&&(I[1] <= border[i].pE[1]+dif)) || ((I[1] <= border[i].pS[1]+dif)&&(I[1] >= border[i].pE[1]-dif))) ) 
				{
					/* �������I �������ཻ�߶�intersectLS֮�䣬�����������˵��е��κ�һ������ô�ͱ���������㡣*/
					FindValidVertex ir;

					/* ���ж� I �� intersectLS[j] �����˵��ľ��룬�������̫�������Ժϲ���һ�� */
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
					
					/* �������I �������ཻ�߶�intersectLS֮�䣬�����������˵��е��κ�һ������ô�ͱ���������㡣*/
					if((I != intersectLS[j].pS) && (I != intersectLS[j].pE))
					{
						ir.pV = I; ir.pVF = I; ir.idx = intersectLS[j].idxS; ir.newIdx = intersectLS[j].idxE; ir.isValid = i; // ����ir.isValid��ʾ�ཻ��������������
						intsResult.push_back(ir);
					}else{
						/* ������� I ���ཻ�߶�intersectLS�������˵�֮һ�غϣ���ôҪ�޸�����˵��״̬��*/
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
					
					/* �������ͶӰ�㣬����������Ӧ�Ŀռ�㡣*/
					/* ����ת���ɾ�γ������ϵ */
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
						/* ���� �ռ�ֱ������ʽ */
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

					/* ���ཻ����뵽��Ч���ϵ�ļ��� */
					ir.pV = I; ir.idx = 0; ir.newIdx = 0; ir.isValid = 0; // ����ir.isValid��ʾ�ཻ��������������
					validBorderVertices.push_back(ir);

					ir.pVF = I; ir.pV = I2; ir.idx = 0; ir.newIdx = 0; ir.isValid = 0; // ����ir.isValid��ʾ�ཻ��������������
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
					
					/* �ҵ�������� �ĸ� borderIntersectLS ���� */
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
						int iFlag = -1;			//�����־��ʶls.pE ��ȡ���Ǳ߽�ε� pS �� ���� pE �㡣
						int addValidBorderFlag = 1;
						LineSeg ls;
						ls.pS = I;   ls.idxS = -1;	ls.validS = 1;
						if(( borderIntersectLS[index].validS == 1) && ( borderIntersectLS[index].validE == 0))
						{	iFlag  = 0; ls.pE = borderIntersectLS[index].pS; ls.idxE = borderIntersectLS[index].idxS;	ls.validE = 1;	}	//ls.validE �ж�����߽�ε������˵����Ļ���Ϣ�ཻ�����2:���ཻ��1:ֻ��һ���ཻ��0:�����ཻ
						if(( borderIntersectLS[index].validS == 1) && ( borderIntersectLS[index].validE == 1))
						{	iFlag  = 0; ls.pE = borderIntersectLS[index].pS; ls.idxE = borderIntersectLS[index].idxS;	ls.validE = 2;	}
						if(( borderIntersectLS[index].validS == 0) && ( borderIntersectLS[index].validE == 1))
						{	iFlag  = 1; ls.pE = borderIntersectLS[index].pE; ls.idxE = borderIntersectLS[index].idxE;	ls.validE = 1;	}
						if(( borderIntersectLS[index].validS == 0) && ( borderIntersectLS[index].validE == 0))
						{	iFlag  = 1; ls.pE = borderIntersectLS[index].pE; ls.idxE = borderIntersectLS[index].idxE;	ls.validE = 0;	}
						ls.border = index;		//border��ʾ��ǰ�ཻ���ǵڼ���

						/* �Ӱ��ཻ�߽�μ����У��Ƴ��뱾���ཻ�����ͬ�Ķ� */
						if(ls.validE == 1)
						{
							for(std::vector<LineSeg>::iterator r = borderHalfIntersectLS.begin(); r != borderHalfIntersectLS.end(); r++)
							{
								if((r->idxS == ls.idxE) && (iFlag == 0)) r->border = -2;
								if((r->idxE == ls.idxE) && (iFlag == 1)) r->border = -2;
							}
						}
						
						/* �����������ͱ߽���ϵ������˵�֮��ľ��룬�������<1.0������Ϊ�ý���Ͷ˵���һ����*/
						/* ע�⣺����ֻ�ж�xy,���ж�z��Ϊ�˺����ǻ��㷨����һ�¡�Galen:2013-02-06 */
						double distS = LenOf2VF(I, borderIntersectLS[index].pS);
						double distE = LenOf2VF(I, borderIntersectLS[index].pE);
						
						/* �������������ڱߵ������˵㣬���Ǻ��Ļ���Ϣ���ཻ�ĵ� ����ls.validE = 2��������������㲻������һ���˵��غϣ� */
						/* �Ǿ�Ҫ�ж������������ߣ���������һ�ߺ��Ļ���Ϣ���ཻ�� */
						/* Ŀǰ�򵥵ز��ò�����ཻ�жϷ���*/
						if((ls.validE == 2) && (distS >= MINDISTTOMERGEVERTEX) && (distE >= MINDISTTOMERGEVERTEX))
						{
							/* ���㽻�����һ���˵�֮��Ĳ�����Ƿ����Ļ���Ϣ�ཻ��*/
							osg::Vec3 IntVer = getIntersectVertexLessTenMeter(I, borderIntersectLS[index].pE);
							osg::Vec3 isect;
							if( TerrainProcessingUtils::getIntersect( gCulture, IntVer, isect ) != false )
							{
								ls.pE = borderIntersectLS[index].pE;	ls.idxE = borderIntersectLS[index].idxE;
							}
							removeValidSegmentWhenWithMultiIntersections.push_back(index);
						}
						
						/* �����ǰ������S�ߵĶ˵�ӽ�һ�¡����ҵ�ǰ�߽�����˵�״̬��ͬ�� */
						if(distS < MINDISTTOMERGEVERTEX) 
						{
							/* ���������ʾ��������Ч�ཻ�� validBorderSegment */
							addValidBorderFlag = 0;

							/* �滻faceVertice �е���ص㡣*/
							ModifyFaceVertexWhenMerge(ir);
							
							/* �ҵ� checkBorderVertices ����ͬ�ĵ㲢���� */
							FindVertexInCBVsets(I);
							
							borderIntersectLS[index].pS = I;						//��Ȼ���Ե���һ���㣬��ô�Ͱѱ߽�εĶ˵㻻������ཻ�㡣
							if (borderIntersectLS[index].validS == 0)	borderIntersectLS[index].validS = 1;
							int idx = index - 1; if(idx <0) idx += borderIntersectLS.size();
							borderIntersectLS[idx].pE = I;				//��Ȼ���Ե���һ���㣬��ô�Ͱѱ߽�εĶ˵㻻������ཻ�㡣
							if (borderIntersectLS[idx].validE == 0)	borderIntersectLS[idx].validE = 1;
							
							/* �жϱ߽綥�㼯terrainBorderVertices �������λ�õ�ԭ���ཻ�ж�ֵnewIdx�Ƿ�Ϊ1 */
							int idxTBV = borderIntersectLS[index].idxS;		//���ֵ��terrainBorderVertices�����Ӧ�±ꡣ
							if(terrainBorderVertices[idxTBV].newIdx >= 0)
							{
								/* �жϵ�ǰindex��borderIntersectLS�߽��ཻ�Σ���index-1�ı߽��ཻ�Σ�����������ཻ�ζ���û�����Ļ���Ϣ�ཻ�ཻ��*/
								/* ˵����ǰ�Ļ���Ϣ���� index �ཻ���ڵ�����߽�����ཻ�����ʱ��Ҫ����validBorderSegment�ཻ�Σ��Ա����ļ��㡣*/
								/* idx = (index - 1)%borderIntersectLS.size() */
								if((borderIntersectLS[idx].validS == 0) && (borderIntersectLS[index].validS == 0) && (borderIntersectLS[index].validE == 0))
								{
									addValidBorderFlag = 1;
								}

								/* �ж�����±��ǰһ��ֵ�Ƿ���Ļ���Ϣ�ཻ��newIdx == 1����*/
								/* ���ǰ���Ǹ�ֵҲ���ཻ����ô��һ���߶�ԭ���Ǹ����ཻ���߶Σ�����֮�������ϳ���ȫ�ཻ���߶Σ��������ǰ��ཻ���߶Σ�����ҪԤ���ų�����*/
								if((idxTBV != 0) && (terrainBorderVertices[idxTBV].newIdx == 0))
								{
									int idxPrev = idxTBV - 1;
									if(terrainBorderVertices[idxPrev].newIdx == 1)
										removeValidSegmentWhenWithMultiIntersections.push_back(idxPrev);
								}

								/* �������±�ֵ��0����Ҫ������������Ǹ�ֵҲ����һ����ֵ����Ϊ����һ���߻���*/
								if(idxTBV == 0)
								{
									int ii = terrainBorderVertices.size() - 1;

									int idxPrev = ii - 1;
									if((terrainBorderVertices[idxPrev].newIdx == 1) && (terrainBorderVertices[ii].newIdx == 0))
										removeValidSegmentWhenWithMultiIntersections.push_back(idxPrev);

									/* ����ֵ�����ֵ���ұ�������϶����Ļ���Ϣ�ཻ����Ϊ�����Ҳ���Ļ���Ϣ�߽�ĵ㣩��*/
									terrainBorderVertices[ii].pVF = I;
									terrainBorderVertices[ii].newIdx = 1;
								}

								/* ����ֵ�����ֵ���ұ�������϶����Ļ���Ϣ�ཻ����Ϊ�����Ҳ���Ļ���Ϣ�߽�ĵ㣩��*/
								terrainBorderVertices[idxTBV].pVF = I;
								terrainBorderVertices[idxTBV].newIdx = 1;		// 1���������ͬ�Ļ���Ϣ�ཻ��0�����������ͬ�Ļ���Ϣ���ཻ��
							}
							
							/* �ཻ���ǽǵ�Ĵ��� */
							CheckSpecialIntVertex(I);
						
							/* ��Ҫ������ǰ���õ�validBorderSegment������û�к�I�غϵĵ㣬����У����滻���� */
							for(unsigned int m = 0; m<validBorderSegment.size(); m++)
							{
								double distVS = LenOf2V(I, validBorderSegment[m].pS);
								double distVE = LenOf2V(I, validBorderSegment[m].pE);
								if(distVS < MINDISTTOMERGEVERTEX) validBorderSegment[m].pS = I;
								if(distVE < MINDISTTOMERGEVERTEX) validBorderSegment[m].pE = I;
							}
							
							/* ����ϲ������� */
							int iFlag = 0;
							for(unsigned int k = 0;k<mergeVertex.size();k++)
							{
								if(index == mergeVertex[k]) {	iFlag = 1; break; }
							}
							if(!iFlag)	mergeVertex.push_back(index);
						}

						/* �����ǰ������E�ߵĶ˵�ӽ�һ�¡����ҵ�ǰ�߽�����˵�״̬��ͬ�� */
						if(distE < MINDISTTOMERGEVERTEX) 
						{
							/* �滻faceVertice �е���ص㡣*/
							ModifyFaceVertexWhenMerge(ir);
							
							/* �ҵ� checkBorderVertices ����ͬ�ĵ㲢���� */
							FindVertexInCBVsets(I);

							borderIntersectLS[index].pE = I;						//��Ȼ���Ե���һ���㣬��ô�Ͱѱ߽�εĶ˵㻻������ཻ�㡣
							if (borderIntersectLS[index].validE == 0)	borderIntersectLS[index].validE = 1;
							unsigned int idx = index + 1; if(idx >= borderIntersectLS.size()) idx = 1;
							borderIntersectLS[idx].pS = I;				//��Ȼ���Ե���һ���㣬��ô�Ͱѱ߽�εĶ˵㻻������ཻ�㡣
							if (borderIntersectLS[idx].validS == 0)	borderIntersectLS[idx].validS = 1;

							/* �жϱ߽綥�㼯terrainBorderVertices �������λ�õ�ԭ���ཻ�ж�ֵnewIdx�Ƿ�Ϊ1 */
							int idxTBV = borderIntersectLS[index].idxE;		//���ֵ��terrainBorderVertices�����Ӧ�±ꡣ
							if(terrainBorderVertices[idxTBV].newIdx >= 0)
							{
								/* ����ֵ�����ֵ���ұ�������϶����Ļ���Ϣ�ཻ����Ϊ�����Ҳ���Ļ���Ϣ�߽�ĵ㣩��*/
								terrainBorderVertices[idxTBV].pVF = I;
								terrainBorderVertices[idxTBV].newIdx = 1;

								/* ���ﻹҪ������β�߽�����⡣*/
								/* ��ΪterrainBorderVertices���һ����͵�һ������һ���㡣*/
								if(idxTBV == terrainBorderVertices.size() -1 )
								{
									terrainBorderVertices[0].pVF = I;
									terrainBorderVertices[0].newIdx = 1;
								}
							}
							addValidBorderFlag = 0;

							/* �ཻ���ǽǵ�Ĵ��� */
							CheckSpecialIntVertex(I);
						
							/* ��Ҫ������ǰ���õ�validBorderSegment������û�к�I�غϵĵ㣬����У����滻���� */
							for(unsigned int m = 0; m<validBorderSegment.size(); m++)
							{
								double distVS = LenOf2V(I, validBorderSegment[m].pS);
								double distVE = LenOf2V(I, validBorderSegment[m].pE);
								if(distVS < MINDISTTOMERGEVERTEX) validBorderSegment[m].pS = I;
								if(distVE < MINDISTTOMERGEVERTEX) validBorderSegment[m].pE = I;
							}

							/* ����ϲ������� */
							int iFlag = 0;
							for(unsigned int k = 0;k<mergeVertex.size();k++)
							{
								if((index+1) == mergeVertex[k]) {	iFlag = 1; break; }
							}
							if(!iFlag)	mergeVertex.push_back(index+1);
						}
						
						if(addValidBorderFlag)	validBorderSegment.push_back(ls);
					}

					/* ���ñ�� */
					bInter = true;
				}
			}
			dif++;
		}
	}//for(j)
	
	/* 	�ж�����һ������������Ч�ཻ����Ϊ0���ཻ��ֻ��һ��������������ǽǵ㣬
		��ô�����ж�����Ļ���Ϣ�����ؾ����ཻ��*/
	if(DetectIntersectStatusOne())
	{
		multiDc.clear();
		return;
	}

	/* ���������� */
	CheckSpecialStatus();

	/* ����Ҫ����һ�����������Ļ���Ϣ������һ�α߽��ж�ε��ཻ������ô�죿*/
	CheckSegmentsIntersection();

	/* �ѱ�������������Ч��Ҳ���뵽validBorderSegment���� */
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

				/* �������߶�֮�䣬�Ƿ��Ѿ�����һ���ཻ�߽��߶�*/
				/* ��鷽���ǣ��ж� validBorderSegment ����� border�������ǰ�κ�border��һ�£�˵���Ѿ�����һ���ཻ�߽���ˡ�*/
				int jFlag = 0;
				for( std::vector<LineSeg>::iterator r= tmpSegment.begin(); r != tmpSegment.end(); r++ )
				{
					if(r->border == borderIdx) { jFlag = 1; break;	}
				}
				if(jFlag) goto nextBorder;
				
				/* ��������Ч�Σ��Ƿ��м��ж�����㣬����еĻ��������������Ч�Ρ� */
				int iFlag = 0;
				for(unsigned int i=0; i<removeValidSegmentWhenWithMultiIntersections.size(); i++)
				{
					if(ls.idxS == removeValidSegmentWhenWithMultiIntersections[i]) { iFlag = 1; break; }
				}
				if(!iFlag)	validBorderSegment.push_back(ls);
			}
nextBorder:
			/* ָ����һ�� */
			borderIdx++;
		}
	}
	
	/* ����Щ�ཻ����pm�߻���ͷ */
	/* Step1:�Ȱ�intsResult ����idx��Ŵ�С����˳������ */
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

	/* Step2:��������intsResult�����뵽fPerimeter��������  */
	for(i=intsResult.size()-1; i>=0; i--)
	{
		j = fPerimeter.size();
		FindValidVertex tmp;
		fPerimeter.push_back(tmp);
		while(j > intsResult[i].newIdx) { fPerimeter[j].pV = fPerimeter[j-1].pV;  fPerimeter[j].isValid = fPerimeter[j-1].isValid; fPerimeter[j].idx = fPerimeter[j-1].idx; fPerimeter[j].newIdx = j; j--;}
		fPerimeter[intsResult[i].newIdx].pV = intsResult[i].pV;
		fPerimeter[intsResult[i].newIdx].isValid = intsResult[i].isValid + 2;  // ��2����ʾ�õ����²�����ཻ�㣬����ʾ�ཻ�����������ϣ�2��ʾ��0��3��ʾ1��...
		fPerimeter[intsResult[i].newIdx].idx = -1;fPerimeter[intsResult[i].newIdx].newIdx = intsResult[i].newIdx;
	}

	/* ��fPerimeter�߻��������Ч�߶�Ҳ����ӵ� validBorderSegment ������ */
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

	/* ���û����Ч�߶Σ���ôֱ���˳�ȥ��*/
	if(validBorderSegment.size() == 0)
	{
		pm->clear();
		newpm->clear();
		return;
	}


	/* �� validBorderSegment ���һ���߻� */
	osg::ref_ptr<osg::Vec3Array> pmLL = new osg::Vec3Array;
	CreateLineLoopFromValidBorderSegment(pmLL);
	
	/* ��� validBorderSegment û���γ���Ч�߻���*/
	if(multiDc.size() == 0)	return;
	if(multiDc[0]->getNumElements() <= 2)
	{
		multiDc.clear();
		return;
	}

	/* Add by Galen 2013-02-18 */
	/* �����ǰ�߻����Ļ���Ϣ�Ǻ�����ô������߻��ϸ�������ӳ��һ�顣���浽 totalVertices ��������*/
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


	/* ���ҵ����ཻ�㣬ȫ�����뵽��Ч����жϼ������档 */
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


/* ���va,vb,vc��������ɵ������Σ��Ƿ�����η���terrain�ཻ��
   ��������дһ���㷨��
   ���ȼ���ؾ���İ�Χ��
   �ټ���ÿ��������İ�Χ��
   ֻҪ��������֮��ľ����������뾶֮�ͣ���˵���治�ཻ������0��
   �����������֮��ľ���С������뾶֮�ͣ���ô���ò�ֵ�İ취����ϸ�жϡ�
   ���������жϺ����õ������㷨�� handleOverlaps()����������ؾ���LINE_LOOP��������LINE_LOOP�Ĺ������������
   �Ͳ���Ҫ������ֵ�㷨�ˡ�
 */
int checkTriwithTerrain(osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vc, osg::Geode *terrain, osg::ref_ptr<osg::Vec3Array> tbLoop)
{
	const osg::BoundingSphere bsTerrain = terrain->getBound();
	osg::ref_ptr<osg::Geode> tri = new osg::Geode;

	/* ���ȣ���Χ����� */	
	/* ���������㶨���һ��LINE_LOOP��Geode���ټ������Χ�� */
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

	/* ����������Χ������֮��ľ��룬��������Χ��뾶֮�� */	
	float sumRadius = bsTerrain.radius() + bsTriangle.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsTriangle.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return 0;
	

	/* ��Σ�������� */
	/* ���������������У����������һ������ͬ�������ཻ���򷵻�1��*/
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

	/* �ٴΣ��������󽻹��� */	
	/* ��������ε�����������ؾ��涼û�н��㣬��ô�п�����������ؾ����н��㡣����������Ǽ�����������ؾ��潻��ĸ�����*/
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

/* ���ȹ�ߵ��������Ƿ�����Ч�㣬ֻҪ��һ��������Ч�㣬�ͷ���0��*/
int checkTriOnSkirt(osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vc)
{     
	/* ������������Ƿ�����Ч�㣬ֻҪ��һ��������Ч�㣬�ͷ���0��*/
	for( std::vector<osg::Vec3>::iterator p= invalidVertices.begin(); p != invalidVertices.end(); p++ )
	{
		if (va==*p || vb==*p || vc==*p)
			return 0;
	}
	return 1;
}


/* ������ϵ������㡣*/
/* �����������, �����������Ļ���Ϣ��, �����������㲻�ڽ�, ��ô����������ǷǷ��ġ� */
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

/* ������ϱ߽���������Ƿ�����Ч�㣬ֻҪ��һ��������Ч�㣬�ͷ���0��*/
int checkTriOnBorder(FindValidVertex *a, FindValidVertex *b, FindValidVertex *c)
{     
	osg::Vec3 va = a->pVF, vb = b->pVF, vc = c->pVF;
	osg::Vec3 sa = a->pV, sb = b->pV, sc = c->pV;
	
	/* ������������Ƿ�����Ч�㣬ֻҪ��һ��������Ч�㣬�ͷ���0��*/
	for( std::vector<osg::Vec3>::iterator p= invalidVertices.begin(); p != invalidVertices.end(); p++ )
	{
		if (va==*p || vb==*p || vc==*p)
		{
			int i=0;
			return 0;
		}
	}

	/* ���ų�������������������ͬ */
	if(( LenOf2V(va,vb) < MINDISTTOMERGEVERTEX ) || ( LenOf2V(va,vc) < MINDISTTOMERGEVERTEX ) || ( LenOf2V(vc,vb) < MINDISTTOMERGEVERTEX ))
	{
		int i=0;
		return 0;
	}

	/* ������������Ƿ��Ǳ߽��ϵĵ㣬����Ǳ߽�㣬��ô��ȡ�ö�Ӧ���������ꡣ*/
	osg::Vec2 ta, tb, tc;
	ta.set(-1.0, -1.0);	tb.set(-1.0, -1.0);	tc.set(-1.0, -1.0);	
	int iBorders = 0;
	for( std::vector<FindValidVertex>::iterator p = allBorderVertices.begin(); p != allBorderVertices.end(); p++ )
	{
		if(va==p->pVF) {	if((ta.x()==-1.0) && (ta.y() == -1.0))	{	ta = p->pT;	iBorders++;	}	};
		if(vb==p->pVF) {	if((tb.x()==-1.0) && (tb.y() == -1.0))	{	tb = p->pT;	iBorders++;	}	};
		if(vc==p->pVF) {	if((tc.x()==-1.0) && (tc.y() == -1.0))	{	tc = p->pT;	iBorders++;	}	};
	}

	/* ����߶ȼ��鷨�� */
	/* ���������������ε���� */
	double area = getArea(sa, sb, sc);
	/* ������һ�� */
	double maxLen = distanceV(sa, sb);
	double lenBC = distanceV(sb, sc);	if(maxLen < lenBC) maxLen = lenBC;
	double lenAC = distanceV(sa, sc);	if(maxLen < lenAC) maxLen = lenAC;
	double heigh = 2 * area / maxLen;
	/* �����������εĸ�̫С�������������������������ڵ��εı��ϡ�����������ξͲ�Ҫ�ˡ�*/
	if( heigh < 1.0 )
	{
		if(iBorders > 1)
			return 0;
		if(( heigh < 0.5 ) && ( iBorders > 0))
			return 0;
	}

	/* �������һ���㲻�Ǳ߽��ϵĵ㣬�򷵻� 1��*/
	if( (ta.x() == -1.0) || (tb.x() == -1.0) || (tc.x() == -1.0) )	return 1;

	/* ��������㶼�Ǳ߽��ϵĵ㣬������������һλ��ͬ��˵������������ͬһ��ı��ϣ��򷵻�0��*/
	if( ((ta.x() == tb.x()) && (ta.x() == tc.x())) || ((ta.y() == tb.y()) && (ta.y() == tc.y())) ) 
	{
		int i=0;
		return 0;
	}

	return 1;
}

#ifdef ADD_SKIRT
/* �����������е�ȹ���� */
std::vector<FindValidVertex> cultSkirtVertice;
std::vector< osg::ref_ptr<osg::DrawElementsUInt> > cultSkirtDrawElements;

void makeSkirt(std::vector<osg::ref_ptr<osg::Vec3Array>> &perimeterList, int currSkirtVertexIndex)
{
	Geodesy geo;
	Carton  cat;
	
	/* ����ȹ��*/
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
			
			/* ���жϳ������,��face���ϵĵ㣨isValid = 1��������culture�ĵ㣨isValid = 0����*/
			for( std::vector<FindValidVertex>::iterator pp= faceVertice.begin(); pp != faceVertice.end(); pp++ )
			{
				fv.isValid = 0;
				if((pp->pVF == fv.pVF) || (pp->pV == fv.pV))
				{
					fv.isValid = 1;
					break;
				}
			}

			/* �ҳ����кϷ�PerimeterVec������� */
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

			/* fv.newIdx ��ʾ����Ч�㼯������û���ҵ�����߻��ϵĵ㡣����ܺ�������Ҳ����֡�*/
			//assert(fv.newIdx != -1);
			if(fv.newIdx != -1)
				fPerimeterlist.push_back(fv);
		}
		if (fPerimeterlist.empty()) continue;
			
		/* ������飬���������ǰһ����ͺ�һ���㶼��culture�ϵĵ㣬���������ȷ��face�ϵĵ㣬��ô���������ĸ̡߳�*/
		unsigned int fsize = fPerimeterlist.size();
		for(unsigned int i=0; i<fsize; i++ )
		{
			/* ��ǰ�� */
			FindValidVertex *fvCurr; fvCurr = &fPerimeterlist[(i+fsize)%fsize];
			/*ǰһ����*/
			FindValidVertex *fvPrec; fvPrec = &fPerimeterlist[(i-1+fsize)%fsize];
			/*��һ����*/
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
				
				/* ���� validVertices ����Ķ�Ӧֵ */
				validVertices[fvCurr->newIdx].pV = fvCurr->pV;
			}
		}
		
		/* �����㼯�Ϻ�DrawElements */
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
   ���㣺 va, vb ������ vO �ľ���֮�͡�
 */
double getDistanceFromTwoVertexToVertex( osg::Vec3 va, osg::Vec3 vb, osg::Vec3 vO)
{

	/* ����va��vO��֮��ľ��롣*/	
	osg::Vec3 Vdisa = getProjectionVertex(vO) - getProjectionVertex(va);
	double Disa = Vdisa.length();

	/* ����vb��vO��֮��ľ��롣*/	
	osg::Vec3 Vdisb = getProjectionVertex(vO) - getProjectionVertex(vb);
	double Disb = Vdisb.length();

	return Disa + Disb;	
}

/* Add by Galen 2013-02-21 
   ��һ��Geometry��ת��Ϊ�ߵļ���
   ���룺osg::Geometry *
   �����std::vector<Edge>
 */
void getEdgesOfGeometry(std::vector<TerrainProcessingUtils::Edge> *kept, osg::Geometry *g)
{
	std::vector<TerrainProcessingUtils::Edge>   edges;
	osg::TriangleFunctor<TerrainProcessingUtils::MakeEdgeList> mel;
	mel.setEdgeVector(&edges);
	g->accept( mel );
	
	/* �ж�edges�������еıߣ��������߳��ֶ�Σ��Ͱ�����߱�עΪfalse����ʾ����������ڲ��ıߡ�2012-10-10 */
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

	/* ��edges�����γ��ֵı��ų�����ֻ����ֻ����һ�εıߣ������һ������ı߽��ϵıߡ�2012-10-10 */
	for( std::vector<TerrainProcessingUtils::Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		if(( p->keep ) && (p->A != p->B ))
		kept->push_back( *p );
	}
}	

 
/* ����ʵ������һ���㷨������һ���������ϸ�����ͬ�Ļ���Ϣ��ı߽�ľ��룬Ҫ��
   1������������ϵĵ����Ļ���Ϣ�治�ཻ��
   2������������ϵĵ�����Ļ���Ϣ������ı��ߵľ���С��һ������ֵ��
 */
void getVextexListNearToCulture(osg::Geode* gdT, osg::Geometry *gmC)
/* ���룺gdT, ӳ���ĵؾ��棬gmC, ӳ�����Ļ���Ϣ�� */
/* ������ؾ����Ϸ��������ĵ㼯�� */
{
	/* Add by Galen 2013-02-22 */
	/* �����жϼ�����g����İ�Χ���Ƿ������gd��Χ�����غ� */
	const osg::BoundingSphere bsTerrain  = gdT->getBound();
	const osg::BoundingSphere bsLakeGmtry = gmC->getBound();

	/* ����������Χ������֮��ľ��룬��������Χ��뾶֮�� */	
	float sumRadius = bsTerrain.radius() + bsLakeGmtry.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsLakeGmtry.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return ;
	/* Add by Galen 2013-02-22 End */
	
	
	/* ����һ���Ļ���Ϣ��Geode, ���������ཻ��*/
	osg::ref_ptr<osg::Geode> flatFeatureGd = new osg::Geode;
	flatFeatureGd->addDrawable(gmC);
	
	/* �����������ת��Ϊ�߻���Ȼ���߻�ת��Ϊ�߼��ϡ�*/
	/* Add by Galen 2013-02-22 ���ж�������ļ������Ƿ��Ѿ��߶λ���*/
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
	
	
	
	/* ���������ÿ�����㡣 */
	for( std::vector<FindValidVertex>::iterator v = faceVertice.begin(); v != faceVertice.end(); ++v )
	{
		double minDistance = 999999.999;
		std::vector<osg::Vec3> minDistVxs;
		
		/* ���ڵ������ϵ�ÿ���㣬����Ҫ����gmC�ཻ */
		osg::Vec3 isect;
		if( TerrainProcessingUtils::getIntersect( flatFeatureGd, v->pVF, isect ))	
		{
			/* Add by Galen 2013-02-25 �������ĸ߳�ֱ�Ӹ������̣߳����ܺ����Ƿ����ù���*/
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

		/* ����߼���kept �������еıߵ������˵㵽 *v ��ľ���֮�͵���Сֵ�����ֵ��С˵�������߾���õ������*/
		for( std::vector<TerrainProcessingUtils::Edge>::iterator p = allLakeEdges.begin(); p != allLakeEdges.end(); p++ )
		{
			double dist = getDistanceFromTwoVertexToVertex(p->A, p->B, v->pVF);
			if(dist < minDistance)
			{
				minDistance = dist;	minDistVxs.clear();	minDistVxs.push_back(p->A);	minDistVxs.push_back(p->B);
			}
		}
		
		/* �������� minDistVxs�� ��� *v ֮��ľ��롣*/
		/* ����߶ȼ��鷨�� */
		/* ���������������ε���� */
		double area = getArea(minDistVxs[0], minDistVxs[1], v->pVF);
		/* ����Ļ���Ϣ�ߵĳ��� */
		double maxLen = distanceV(minDistVxs[0], minDistVxs[1]);
		double heigh = 2 * area / maxLen;
		
		/* ��¼��������Ϣ */
		if(heigh < MINDISTANCETOEDGE)
		{
			VertexDistance vd;	vd.pVF = v->pVF;	vd.distance = heigh;	vertexNearLake.push_back(vd);
		}
	}	// v

}
/* Add by Galen 2013-02-20  End */


/* �ع�features, �� features �ϲ����ؾ��ཻ����ɾ����
   Ŀǰ����㷨������̫���ˡ�
   ��������дһ���㷨��
   ���ȼ���ؾ���İ�Χ��
   �ټ���ÿ��������İ�Χ��
   ֻҪ��������֮��ľ����������뾶֮�ͣ���˵���治�ཻ��
   �����������֮��ľ���С������뾶֮�ͣ���ô���ò�ֵ�İ취����ϸ�жϡ�
 */
bool reConstructFeature(osg::Geode* gd, osg::Geometry *g, osg::ref_ptr<osg::Vec3Array> tbLoop)
{
	/* �����жϼ�����g����İ�Χ���Ƿ������gd��Χ�����غ� */
	const osg::BoundingSphere bsTerrain  = gd->getBound();
	const osg::BoundingSphere bsLineloop = g->getBound();

	/* ����������Χ������֮��ľ��룬��������Χ��뾶֮�� */	
	float sumRadius = bsTerrain.radius() + bsLineloop.radius();
	osg::Vec3 Vdis = bsTerrain.center() - bsLineloop.center();
	float centerDis = Vdis.length();
	if(centerDis > sumRadius) return false;

	/* ���жϼ�����g��ÿ�������εİ�Χ���Ƿ������gd��Χ�����غ�*/
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

/* ����һ�������(�߻�)�������*/
double GetAreaOfLineloop(osg::ref_ptr<osg::Vec3Array> ll)
{
	/* �Ȱѵ�ǰperlist�߻������ǻ���һ���棬Ϊ�˸��߽�μ������ཻ��� */				
	osg::ref_ptr<osg::Geometry> g = ConvertFromLINELOOPToTranglesFace(ll);
	if(g == NULL)	return 0.0;
	
	/* �������е�����������֮�͡�*/
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


/* ��һ�������ж��primetives��dc�ϲ���dclist���� */
void MergeDcIntoDclist(std::vector<osg::ref_ptr<ArealConstraint> > *dclist, osg::ref_ptr<ArealConstraint> mergedc, string dcname)
{
	/* �Ȱ����ϲ���һ��gmlist */
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

	/* ����Ǻ�������ô����ÿƬ��������������ĳ�����������С��С���������ɺ��Բ��� */
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
	
	/* ���߻����dcѹ�뵽dclist */
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


/* �����еľ�����ͬ���Ƶ�dc�ϲ���Ϊһ���������㸲����� */
void MergeDclistStepOne(std::vector<osg::ref_ptr<ArealConstraint> > *dclist)
{
	if(dclist->size() == 0)	return;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/* ������δ������ �Ļ���Ϣ·�н��������Ĵ���*/
	/* ��ͳ���ж��ٸ�ͬ���� dc������1�����ϵģ����� merge */
	totalRoad = 0;
	totalIsland = 0;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist->begin(); p != dclist->end(); p++ )
	{
		if (p->get()->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING) totalRoad++;
		if (p->get()->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING) totalIsland++;
	}

	/* �����е� ����Ϊ Road �ĺϲ���һ�� Road */
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
	
		/* ��mergedc�����ҳ������Ķ��㡣*/
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
			
			/* �ϲ�mergedc */
			if(mergedc.valid())		MergeDcIntoDclist(dclist, mergedc, "Road");
			if(mergedc2.valid())	MergeDcIntoDclist(dclist, mergedc2, "Road");
		}
	}	//totalRoad > 1;
}


/* ��Lake��Island��dc�ϲ���Ϊһ���������㸲����� */
void MergeLakeAndIsland(std::vector<osg::ref_ptr<ArealConstraint> > *dclist)
{
	if(dcLakes.size() == 0) return;
	
	/* ��ͳ���ж��ٸ�"Island" ���� ��Lake�� dc������1�����ϵģ����� merge */
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
	/* ������δ������ �Ļ���Ϣ ���������ɵ��� ������Ĵ���*/
	/* ��Ҫ�ڵ��Ļ���Ϣriver�� ��Ҫ�ڵ� island ֮��Ҳ�� merge*/
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



/* ����ͬ��̬�ݵ����ཻ�ĵ������ */
void CreateValidVertexOnFaceOfDynamicGrass(MaxMinLonLat *mmll, osg::ref_ptr<osg::Geometry> gm, double elevMid, double elevOfGrass)
{
	/* ����� 1 ����ռ�ݵľ�γ���Ƕ��� deltaLon, deltaLat */
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

	/* ������̬�ݽ�� */
	osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
	gnode->setName("DynamicGrass");
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry;
	gnode->addDrawable(geom);
	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	geom->setVertexArray(vx);

	/* ����� */
	for(double i = mmll->minLong; i < mmll->maxLong; i += UNITDISTANCE * deltaLon)
	{
		for(double j = mmll->minLati; j < mmll->maxLati; j += UNITDISTANCE * deltaLat)
		{
			/* ��γ��ת����xyz */
			geo._lon = i; geo._lat = j; geo._elevation = elevMid;	GeodToCart(&geo, &car);
			osg::Vec3 currVert(car.x, car.y, car.z);

			/* ����ÿ�����Ƿ�� g ���ཻ���󽻵� */
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

/* ������̬�ݵ��� */
void CreateFaceOfDynamicGrass(osgTDS::TexCoordPreserver *tcp)
{
	if(dclist.empty()) return;

	//����һ���µĵ㼯�ϣ�ֻ�������ϵĵ�
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
		
		/* �ҳ���ǰ��̬�ݵĸ߶� */
		SplitName(name);
		if(NameParts.size() >= 2)	elevOfGrass = atof(NameParts[1].c_str());

		if (p->get()->getinteriorTrisSize()>0)
		{
			if(1)
			{
				osg::ref_ptr<osg::Vec3Array> arpts=p->get()->getPoints(dcoords.get());

				/* �����ǻ�������е㼯�϶�����validVertices���� */
				/* ����validVertices�������е㣬��ͶӰ��ӳ�䵽�ռ�� */
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

				//�ҳ����е���Ч�档����� validElement ����
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

				//����һ���µĵ����У�����������Ч��
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//��ֻ�����Ϸ������ε�validElement�����Ԫ�ؼӽ���
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

				/* ������������и��ļ����档��ȫ�ռ�������档*/
				osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
				gm->setVertexArray( validCoords.get() );
				gm->addPrimitiveSet( dui );;

				/* �ҳ��������и�󼸺����ཻ�ĵ㡣�����ļ����2��*/
				/* ���ݼ�������Ķ�̬�ݵ��棬������ڲ��Ĳ���� */
				CreateValidVertexOnFaceOfDynamicGrass(&mmll, gm, elevMid, elevOfGrass);
				
			}
		}
	} 
}



/* ��ϸ�������dc�����ó����� */
void MoveOutDetailTxtDc()
{
	/* �ȰѺ����ó�����*/
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
	
	
	/* �ٰ�ϸ���������ó��� */	
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


/* ���������� */
void CreateFaceOfLake(osgTDS::TexCoordPreserver *tcp)
{
	if(dcLakes.empty()) return;


	//����һ���µĵ㼯�ϣ�ֻ�������ϵĵ�
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
		
		/* ���name��β�пո����������� */
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
				/* �ҳ���Ч���㼯�� */
				osg::ref_ptr<osg::Vec3Array> arpts=p->get()->getPoints(dcoords.get());

				/* �����ǻ�������е㼯�϶�����validVertices���� */
				/* ����validVertices�������е㣬��ͶӰ��ӳ�䵽�ռ�� */
				int j = 0;
				validVertices.clear();	validElement.clear();
				for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
				{
					FindValidVertex vVert;
					osg::Vec3 pSpace = getSpaceVertex(*r);
					vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
					validVertices.push_back(vVert);
				}

				//�ҳ����е���Ч�档����� validElement ����
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

				//����һ���µĵ����У�����������Ч��
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//��ֻ�����Ϸ������ε�validElement�����Ԫ�ؼӽ���
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

#ifdef ADD_SKIRT
				/* ����ȹ�� */
				CreateSkirt(validCoords, dui, 0);
#endif

				/* ������������и��ļ����塣*/
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

/* ����ϸ��������棬���Ļ���Ϣ�� */
void CreateFaceOfDetailTxtWhenNoCulture(osgTDS::TexCoordPreserver *tcp)
{
	if(dcDetail.empty()) return;

	//����һ���µĵ㼯�ϣ�ֻ�������ϵĵ�
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
		
		/* ���name��β�пո����������� */
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

				/* �����ǻ�������е㼯�϶�����validVertices���� */
				/* ����validVertices�������е㣬��ͶӰ��ӳ�䵽�ռ�� */
				int j = 0;
				validVertices.clear();	validElement.clear();
				for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
				{
					FindValidVertex vVert;
					osg::Vec3 pSpace = getSpaceVertex(*r);
					vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
					validVertices.push_back(vVert);
				}

				//�ҳ����е���Ч�档����� validElement ����
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

				//����һ���µĵ����У�����������Ч��
				osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
				for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
					validCoords->push_back(pp->pV);

				//��ֻ�����Ϸ������ε�validElement�����Ԫ�ؼӽ���
				osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
				for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
					dui->push_back(*p);

				if (dui->empty())	return;

#ifdef ADD_SKIRT
				/* ����ȹ�� */
				CreateSkirt(validCoords, dui, 0);
#endif

				/* ������������и��ļ����塣*/
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

/* ����ϸ��������棬���Ļ���Ϣ�� */
void CreateFaceOfDetailTxtWhenHasCulture(osg::ref_ptr<osg::Geometry> geom)
{
	if(dcDetail.empty()) return;

	osg::ref_ptr<osg::Geode> gnode = NULL;
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
	{
		std::string name= p->get()->getName();
		
		/* ���name��β�пո����������� */
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

/* ���dclist������ǲ��Ǻ����͵��������ͬ��dc����������Ҫ�� 
   1��dclistֻ������dc��
   2������dc����������һ���ġ�
   3������dc������Ҳһ����
   ��������������һ��dc��ֻ��һ����
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
	
	/* �������������ͬ������Ҳ��ͬ */
	if(dc0->getName() == dc1->getName())
	{
		dclist->clear();
		osg::ref_ptr<ArealConstraint> dc = dynamic_cast<ArealConstraint *> (dc0.get());
		dclist->push_back(dc);
	}

	return ret;
}

/* ��ʼ�������õ������� */
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

/* �������ؾ����ϵ����е��Ϊ�������֣�һ�����ǵؾ����ϵĵ㣬һ������ȹ���ϵĵ㡣*/
/* ����������������ַ�������Ϣ��faceVertices��skirtVertices���ַ�ԭ���ǣ�����һ����������Ϣ��ͬ���߶ȵ�����һ����ģ�����skirt�㡣 */
void SplitTerrainVertexToFaceAndSkirt()
{
	/* ���ȵ�һ����϶���face���ϵĵ㡣 */
	faceVertice.push_back(terrainVertices.front());
	
	/* �ٰ�������������ͬ������߶ȸߵ��������� ȹ�ߵ� ���� ���ϵ� */
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

/* ���߽��߻��ϵ����е㣬��ȹ�ߵ㣬���滻�����ϵĵ� */
void ChangeBorderVertexFromSkirtToFace(osg::Geometry *g)
{
	std::vector<osg::ref_ptr<osg::Vec3Array> > periList;

	/* ������ؾ���Ե���߻�����Щ���п��ܶ���ȹ���ϵĵ㡣*/
	TerrainProcessingUtils::computePerimeterAll(g, periList);
	if(periList.size() > 0) terrainBorderLineLoop = periList[0];


	/* ��border�߻��ϵ����ж�����Ϣ���뵽 terrainBorderVertices ���顣*/
	unsigned int j = 0;
	for( osg::Vec3Array::iterator r = terrainBorderLineLoop->begin(); r != terrainBorderLineLoop->end(); r++ )
	{
		FindValidVertex fvv;
		fvv.pV = *r; fvv.idx = j; fvv.newIdx = 0; j++; fvv.isValid = -1;
		terrainBorderVertices.push_back(fvv);
	}
    
	/* ���� terrainBorderVertices ���������pVֵ �� skirtVertics ��pVֵһһ��Ӧ���ҳ� terrainBorderVertices �����Ӧ�� pTֵ ��*/
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
	{
		for( std::vector<FindValidVertex>::iterator f = skirtVertics.begin(); f != skirtVertics.end(); f++ )
		{
			if(f->pV == p->pV) {	p->pT = f->pT; break;	}
		}
	}
	
	/* �ٸ��� terrainBorderVertices ���������pTֵ �� faceVertics ��pTֵһһ��Ӧ���滻���� terrainBorderVertices �����Ӧ�� pV �� pVF ֵ ��*/
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
	{
		for( std::vector<FindValidVertex>::iterator f = faceVertice.begin(); f != faceVertice.end(); f++ )
		{
			if(f->pT == p->pT) {	p->pV = f->pV; p->pVF = f->pVF; break;	}
		}
	}
	
	/* ���� terrainBorderVertices ���������pTֵ �ҳ��ĸ��ǵ��λ�ã���conner���±�һ�£������� isValid����, ���඼��-1 */
	for( std::vector<FindValidVertex>::iterator f = terrainBorderVertices.begin(); f != terrainBorderVertices.end(); f++ )
	{
		if((f->pT[0] == 0.0) && (f->pT[1] == 0.0))  { f->isValid = 0;}
		if((f->pT[0] == 0.0) && (f->pT[1] == 1.0))  { f->isValid = 1;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 1.0))  { f->isValid = 2;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 0.0))  { f->isValid = 3;}
	}
	
	/* �ٻ�ԭ�� terrainBorderLineLoop */
	terrainBorderLineLoop = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator p = terrainBorderVertices.begin(); p != terrainBorderVertices.end(); p++ )
		terrainBorderLineLoop->push_back(p->pV);

	/*  ����߽���ÿһ����ľ�γ��ӳ�䡣 */
	borderPMs.clear();
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
	{
		ProjectionMap pm;
		getProjectionMapByVertex(&pm, q->pV);
		borderPMs.push_back(pm);	
	}
}

/* �Ӵ�ȹ�ߵĸ߶� 2013-01-04 */
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

	/* ���˸߶ȵ�ҲҪ����д�ؾ� */
	CurrentTerrainIsCut = 1;
}



/* ����csmList������faceVertice������ĸ߶ȣ��� coords��������ĸ߶� */
void TuneElevationsByCSMList(osg::ref_ptr<osg::Vec3Array> coords)
{
	if(csmList.size())
	{
		/* �ȼ����Χ���Զ����ϵ*/
		for(std::vector< osg::ref_ptr<CorrectionSpaceMesh> >::iterator c = csmList.begin(); c != csmList.end(); c++)
		{
			osg::ref_ptr<CorrectionSpaceMesh> csm = *c;

			/* ����������Χ������֮��ľ��룬��������Χ��뾶֮�� */	
			float sumRadius = terrainBS.radius() + csm->radius;
			osg::Vec3 Vdis = terrainBS.center() - csm->center;
			float centerDis = Vdis.length();
			if(centerDis < sumRadius) validCSM.push_back(csm);
		}
		if(validCSM.size() == 0) return;
			
#if 1		//�����Χ��İ취�е�ȱ�ݡ����һ�������ĵ��ξͲ��ʺ��ˡ�
		/* Add by Galen 2013-03-15 �����°취������csmList����������Ƿ��ཻ��*/
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

		/* ���˸߶ȵ�ҲҪ����д�ؾ� */
		CurrentTerrainIsCut = 1;
		
		/* �ȵ����������ϵĵ�ĸ߶� */
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

		/* �ٵ�������ȹ���ϵĵ�ĸ߶� */
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

/* ��faceVertice�������ҳ� �ǲ� �� �߲� ��һЩ��Ϣ */
void FindCornerAndBorderInfoInFaceVertices()
{
	/* ��faceVertice�����в����ĸ��ǵĶ�����Ϣ���Ӷ��ó��ĸ��ߵ���Ϣ  */
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
	{
		if((f->pT[0] == 1.0) || (f->pT[0] == 0.0) || (f->pT[1] == 1.0) || (f->pT[1] == 0.0))
		{
			FindBorderVertex fbv;	fbv.fvv.pV = f->pV; fbv.fvv.pVF = f->pVF; fbv.fvv.newIdx = f->newIdx; fbv.fvv.idx = f->idx; fbv.fvv.isValid = 0; fbv.onDetailed = 0; fbv.onCulture = 0;
			checkBorderVertices.push_back(fbv);
			f->isValid = 2; // ������ʾ���ϵĵ���border�ϵĵ�
		}else{
			FindValidVertex fvv;	fvv.pV = f->pV; fvv.pVF = f->pVF; fvv.newIdx = f->newIdx; fvv.idx = f->idx; fvv.isValid = 0;
			terrainInnerVertices.push_back(fvv);
		}

		if((f->pT[0] == 0.0) && (f->pT[1] == 0.0))  { conner[0].pV = f->pV; conner[0].pVF = f->pVF; conner[0].pT = f->pT; conner[0].idx = f->idx; conner[0].newIdx = 0; conner[0].isValid = 0;}
		if((f->pT[0] == 0.0) && (f->pT[1] == 1.0))  { conner[1].pV = f->pV; conner[1].pVF = f->pVF; conner[1].pT = f->pT; conner[1].idx = f->idx; conner[1].newIdx = 1; conner[1].isValid = 0;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 1.0))  { conner[2].pV = f->pV; conner[2].pVF = f->pVF; conner[2].pT = f->pT; conner[2].idx = f->idx; conner[2].newIdx = 2; conner[2].isValid = 0;}
		if((f->pT[0] == 1.0) && (f->pT[1] == 0.0))  { conner[3].pV = f->pV; conner[3].pVF = f->pVF; conner[3].pT = f->pT; conner[3].idx = f->idx; conner[3].newIdx = 3; conner[3].isValid = 0;}
	}

	/* ���㽫�ĸ������ӵ�һ���Ǽ��� */	
	calcConnerEx();

	/* �������ĸ�����ɵļ������� */
	border[0].pS = conner[0].pVF; border[0].pE = conner[1].pVF;
	border[1].pS = conner[1].pVF; border[1].pE = conner[2].pVF;
	border[2].pS = conner[2].pVF; border[2].pE = conner[3].pVF;
	border[3].pS = conner[3].pVF; border[3].pE = conner[0].pVF;

	/* ����һ���򵥱߻������������ع���ʱ��ʹ�õ� */
	terrainBoundLoop = new osg::Vec3Array;
	for (unsigned int ic=0;ic<4;ic++) terrainBoundLoop->push_back(conner[ic].pVF);
		terrainBoundLoop->push_back(conner[0].pVF);
}

/* ���ݵ����Ϣ��ȷ�������Ϣ����������������㶼�����ϵĵ㣬���������ǵؾ��棬�������ȹ���档 */
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

/* ��һ��LINE_LOOP��ת���� ȹ�� ��geometry */
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToSkirtFace(osg::ref_ptr<osg::Vec3Array> gvx)
{
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	
	/* ���ָ����⣬����߻�ֻ���ĸ��㣬Ҳ����ֻ��һ�������ε�ʱ�򣬵��������ǻ�����������֧�������㡣2012-10-10*/
	if(gvx->size() >= 4)
	{
		/* �������߻���һ��һ�������߻���*/
		osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;
		for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
		{
			coords->push_back( *pvx );
			osg::Vec3d gI = cTog( *pvx );	gI.z() -= SKIRT_HEIGHT;
			osg::Vec3 cI = gToc( gI );
			coords->push_back( cI );
		}
		gm->setVertexArray(coords);

		/* �ٲ����� */
		/* ���� PrimitiveSet����ͬ����Ҫ�в�ͬ������˳��*/
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


/* ��һ��LINE_LOOP��ת���� �������� ��geometry */
osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToTranglesFace(osg::ref_ptr<osg::Vec3Array> gvx)
{
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	
	/* ���ָ����⣬����߻�ֻ���ĸ��㣬Ҳ����ֻ��һ�������ε�ʱ�򣬵��������ǻ�����������֧�������㡣2012-10-10*/
	if(gvx->size() > 4)
	{
		/* �������߻���������ͬ���߻�*/
		osg::ref_ptr<osg::Vec3Array> perimeter1 = new osg::Vec3Array;
		osg::ref_ptr<osg::Vec3Array> perimeter2 = new osg::Vec3Array;
		for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
		{
			perimeter1->push_back( *pvx );
			perimeter2->push_back( *pvx );
		}

		/* ������һ���߻������������ǻ�����һ���߻���Ϊ������Լ���� */
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

/* ��һ��dclist�����ҳ�������Ч�㣨���ǲ�������߻��ĵ㣩�ļ��ϣ�*/
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


/* ��һ��dc��ת���� �������� ��geometry */
osg::ref_ptr<osg::Geometry> ConvertFromDCToTranglesFace(osg::ref_ptr<ArealConstraint> dc)
{
	osg::ref_ptr<osg::Vec3Array> vtx = FindValidVertexInDC(dc);
	if( vtx == NULL ) return NULL;

	/* ʹvtx�߻���֤��β������ */
	if (!vtx->empty() && vtx->front()!=vtx->back())
		vtx->push_back(vtx->front());

	osg::ref_ptr<osg::Geometry> gm = ConvertFromLINELOOPToTranglesFace(vtx);
	return gm.release();
}

#ifdef ADD_SKIRT
	/* ����ȹ�� */
void CreateSkirt(osg::ref_ptr<osg::Vec3Array> validCoords, osg::ref_ptr<osg::DrawElementsUInt> dui, int mode)
{
	/* ���ȳ�ʼ��ȫ�ֱ��� */
	cultSkirtVertice.clear();
	cultSkirtDrawElements.clear();
	
	/* ������α߽��н��� �� ����α߽��޽��� ����������ֱ����ȹ�� */
	if (intersectlist.empty())
	{
		if(mode == 1)
		{
			//��ԭ��ȹ�ߵĵ�Ҳ���뵽�㼯�ϵ��У���ȷ��ȹ�ߵ������newIdx
			int startIdx = validCoords->size();
			for( std::vector<FindValidVertex>::iterator p = skirtVertics.begin(); p != skirtVertics.end(); p++ )
			{
				validCoords->push_back(p->pV); 
				p->newIdx = startIdx++;
			}
	
			//ԭ����������֮ǰ�����ϵĵ㼯�ϣ���Ҫ����Ч�㼯��validVertices������ҳ������˵�����֮����µ�idx;
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
	 
			//��ȹ����Ԫ�ؼ�skirtElement�����ϵ�idx�滻���µ�newIdx��Ϣ
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
	
			//����������ȹ��skirt����
			for( std::vector<unsigned int>::iterator p = skirtElement.begin(); p != skirtElement.end(); p++ )
				dui->push_back(*p);
				
			/* �Ȱ�ȹ�ߵĵ�Ҳ���뵽validVertices�������� */
			for( std::vector<FindValidVertex>::iterator p = skirtVertics.begin(); p != skirtVertics.end(); p++ )
			{
				FindValidVertex fvv;
				fvv.pV = p->pV; fvv.pVF = p->pVF; fvv.pT = p->pT; fvv.idx = p->idx; fvv.newIdx = p->newIdx; fvv.isValid = p->isValid;
				validVertices.push_back(fvv);
			}
	
			/* ���ݱ߽�ζ˵����Ч�ԣ�������ȹ���ǲ��ǿ��Ի�������*/
			//�ҳ����е���Ч�档����� validElement ����
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
			//��ֻ�����Ϸ������ε�validElement�����Ԫ�ؼӽ���
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

	/* ����ȹ��֮ǰ�����һ���߻�������Ƿ���ʡ�*/
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

	/* ����ȹ�� */
	makeSkirt(perimeterList, validCoords->size());

	/* ������makeSkirt()�����п����޸���validVertices�����ֵ��������Ҫ���´���validCoords */
	validCoords->clear();
	for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
		validCoords->push_back(pp->pV);

	/* �����������õ�ȹ�ߵ� */
	for( std::vector<FindValidVertex>::iterator p = cultSkirtVertice.begin(); p != cultSkirtVertice.end(); p++ )
		validCoords->push_back(p->pV);
}

#endif	//ADD_SKIRT


/* ����ÿһ��Geometry */
osg::ref_ptr<osg::Geometry> processGeometry( osg::Geometry *g, const GeodeList &features )
{
	elevationOfLake = -9999.99;
	pCurTerrain = g;
	pCurrTerrainGeode = new osg::Geode;
	pCurrTerrainGeode->addDrawable(g);
	osg::ref_ptr<osg::Vec3Array> coords    = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	osg::ref_ptr<osg::Vec2Array> texCoords = dynamic_cast<osg::Vec2Array *>(g->getTexCoordArray(0));
	osg::ref_ptr<osg::Vec3Array> norCoords = dynamic_cast<osg::Vec3Array *>(g->getNormalArray() );

	/* �ؾ�����Ч�Լ�� */
	if( !coords.valid() )    return NULL;
	if( coords->size() == 0) return NULL;
    
	/* ���㵱ǰ�ؾ���İ�Χ�� */
	terrainBS = g->getBound();
    
	/* ��һ������ȹ��skirt�����ԭGeometry���ų���ȥ�� */
	/* �´���һ��Geometry���������������g�����ȹ��skirt����档*/
	/* ��ʼ�������õ������� */
	InitVectors();

	/* �������ж�����Ϣ��������Ϣ��terrainVertices����, �������ͶӰ����ƽ��ĸ������ֵ�� */
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
    
	/* �������ؾ����ϵ����е��Ϊ�������֣�һ�����ǵؾ����ϵĵ㣬һ������ȹ���ϵĵ㡣*/
	SplitTerrainVertexToFaceAndSkirt();

	/* ���û��ȹ�ߣ����ǵؾ����п����� features �ļ��в���ʶ���gnode */
	if (skirtVertics.empty()) return NULL;

#ifdef ADD_SKIRT
	/* ��������ȹ�ߵĸ߶� */
	Geodesy fgeo, sgeo;
	Carton  fcat, scat;
	fcat.x = faceVertice[0].pV[0];	fcat.y = faceVertice[0].pV[1];	fcat.z = faceVertice[0].pV[2];
	scat.x = skirtVertics[0].pV[0];	scat.y = skirtVertics[0].pV[1];	scat.z = skirtVertics[0].pV[2];
	assert(skirtVertics[0].pT==faceVertice[0].pT);
	CartToGeod(&fcat, &fgeo);	CartToGeod(&scat, &sgeo);
	skirtLength = fgeo._elevation - sgeo._elevation;
#endif

	/* ����ȹ�ߵĵ㣬�Ӵ�ȹ�ߵĸ߶� 2013-01-04 */
#ifdef	MAKE_HIGH_SKIRT
	skirtLength *= 4;
	AddSkirtHeight();
#endif

	/* ����csmList������faceVertice��coords��������ĸ߶� */
	TuneElevationsByCSMList(coords);


	/* ��ʼ���������걣������ */
	osgTDS::TexCoordPreserver tcp(&terrainVertices);
    	
	/* ��faceVertice�������ҳ� �ǲ� �� �߲� ��һЩ��Ϣ */
	FindCornerAndBorderInfoInFaceVertices();

	/* ���߽��߻��ϵ����е㣬��ȹ�ߵ㣬���滻�����ϵĵ� */
	ChangeBorderVertexFromSkirtToFace(g);
	
	/* ���ݵ����Ϣ��ȷ�������Ϣ����������������㶼�����ϵĵ㣬���������ǵؾ��棬�������ȹ���档 */
	SplitFaceInfoByVertex(g);
	
	// ���������������ᣬֻ�����ϵĵ㣬�ų�ȹ�ߵĵ㡣
	// �ȴ���һ��Geode��ֻ�������ϵ������Σ��㼯�Ͼ�ȫ������������Щ�㶼ת����ͶӰ�㣬ͶӰ����ƽ�档(Ҳ���ǵ�ת���ɾ�γ������ϵ�󣬸߶�ֵ����Ϊ0.0)�������������ᡣ
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

		/* �Ȱѵ��η����ĸ������㣬Ԥ�����㵱ǰfeatures����η����ĸ��ǵ��ཻ�����*/
		for (ic=0;ic<4;ic++)	conner[ic].isValid = 0;

		osg::Geode* f = (*p)._geode.get();
		
		/* ���涨���flatFeatureGd, ���Ļ���Ϣ�����е㶼����ͶӰ��Geode��*/
		osg::ref_ptr<osg::Geode> flatFeatureGd = new osg::Geode;
		for( unsigned int ii = 0; ii < f->getNumDrawables(); ii++ )
		{
			osg::ref_ptr<osg::Geometry> g = dynamic_cast<osg::Geometry *>(f->getDrawable(ii));
			osg::ref_ptr<osg::Geometry> fg = new osg::Geometry(*g, osg::CopyOp::DEEP_COPY_ALL);
			std::string fname;
			if (g->getName().empty()==false)  fname = g->getName(); else fname = f->getName();

			
			/* �����ǰ������Lake����ô��ȡ����ĸ߶� */
			if(fname.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			{
				if(elevationOfLake < 0.0)
				{
					osg::ref_ptr<osg::Vec3Array> ver =  dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
					osg::Vec3Array::iterator q = ver->begin();
					
					/* ȡ����һ����ĸ߶ȼ��ɡ�*/
					Geodesy geo;
					Carton  cat;
	
					cat.x = q->x(); cat.y = q->y(); cat.z = q->z();	CartToGeod(&cat, &geo);	
					double currLakeElev = geo._elevation;
					
					if(currLakeElev > elevationOfLake) elevationOfLake = currLakeElev;
				}
			}

			/* �����ǰ������Airport����ô��ȡ�����ĸ߶� */
			isAirport = false;
			m_elevOfAirport = 0.0;
			if((fname.substr(0, strlen(DEFINE_AIRPORT_STRING)) == DEFINE_AIRPORT_STRING) || (fname.substr(0, 7) == "AIRPORT"))
			{
				SplitName(fname);
				if(NameParts.size() >= 2)	m_elevOfAirport = atof(NameParts[1].c_str());
				isAirport = true;
			}

			/* ��ȡ���������*/
			osg::Vec3Array *vert =  dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
			if(vert == NULL)	continue;
			int vert_size =  vert->size();

			/* ����Geometry����ΪLINE_LOOP��Ҫ��LINE_LOOP�ָ���������͡� */
			osg::ref_ptr<osg::PrimitiveSet> ps = g->getPrimitiveSet(0);
			if(ps->getMode() == osg::PrimitiveSet::LINE_LOOP)
			{
				osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
				fg = ConvertFromLINELOOPToTranglesFace(gvx);
				if(fg == NULL)	continue;
			}
			
			/*��fg���ж��㶼��ͶӰ������ͶӰ����ƽ�档*/
			makeGeometryFlat(fg);
			flatFeatureGd->addDrawable(fg);

			/* �ع�features֮ǰ������Ҫ����һ�µؾ��ĸ�����features���ཻ�����
			   ��Ϊ�������˵ؾ�һ���ǣ���features�ཻ����ȴ���ع����features���ཻ ���������ɴ���
			   �ؾ�������ĸ��ǵ���features�Ƿ��н��㣬�н������isVaild = 1  */
			osg::Vec3 isect;
			for (ic=0;ic<4;ic++)
			{
				if( TerrainProcessingUtils::getIntersect( flatFeatureGd, conner[ic].pVF, isect ))	conner[ic].isValid += 1; // ��ʾ�ǵ�� feature �ཻ
			}

			/* �鿴�������еĵ��Ƿ�� feature �ཻ�� */
			int TotalIntersectVerts = 0;
			for( std::vector<FindValidVertex>::iterator q= faceVertice.begin(); q != faceVertice.end(); q++ )
			{
				if (TerrainProcessingUtils::getIntersect( flatFeatureGd, q->pVF, isect ))	
				{	
					subIntersectCt++;	
					TotalIntersectVerts++;	
					
					/* ����Ļ���Ϣ���ǻ����棬*/
					if((isAirport == true) && (m_elevOfAirport != 0.0))
					{
						Geodesy geo;
						Carton  cat;
					
						cat.x = q->pVF.x(); cat.y = q->pVF.y(); cat.z = q->pVF.z();
						CartToGeod(&cat, &geo);	geo._elevation = m_elevOfAirport; GeodToCart(&geo, &cat);
						q->pV.x() = cat.x; q->pV.y() = cat.y; q->pV.z() = cat.z;
						
						/* ���ø��¸̵߳ı�־���� */
						UpdateTheElevationInTerrain = 1;
					}
				}
			}
			
#ifdef ADD_SKIRT
			/* ����һ�£��ٲ鿴ȹ�������еĵ��Ƿ�� feature �ཻ��*/
			for( std::vector<FindValidVertex>::iterator q= skirtVertics.begin(); q != skirtVertics.end(); q++ )
			{
				if (TerrainProcessingUtils::getIntersect( flatFeatureGd, q->pVF, isect ))	
				{	
					/* ����Ļ���Ϣ���ǻ����棬*/
					if((isAirport == true) && (m_elevOfAirport != 0.0))
					{
						Geodesy geo;
						Carton  cat;
					
						cat.x = q->pVF.x(); cat.y = q->pVF.y(); cat.z = q->pVF.z();
						CartToGeod(&cat, &geo);	geo._elevation = m_elevOfAirport - skirtLength; GeodToCart(&geo, &cat);
						q->pV.x() = cat.x; q->pV.y() = cat.y; q->pV.z() = cat.z;

						/* ���ø��¸̵߳ı�־���� */
						UpdateTheElevationInTerrain = 1;
					}
				}
			}
#endif
			
			/* ����ؾ��������еĵ㶼���Ļ���Ϣ���ཻ����ô���԰��Ļ���Ϣ��򵥴���һ�¡�ֱ�Ӹ�������ؾ���Ϳ��ԡ�*/
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
			/* �����ǰ�Ļ���Ϣ�Ǻ�������ôҪ���ٽ����ߵĶ�����һЩ�ռ����Ա������ƽ������*/
			if(fname.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)
			{
				getVextexListNearToCulture(flatTerrainGd.get(), gg.get());
			}
			/* Add by Galen 2013-02-21 End */

			/* �ع�features, �� features �ϲ����ؾ��ཻ����ɾ��������·���ؾ��ж���ཻ����� */
			/* ������ͶӰ��ĵؾ��� �� ͶӰ����Ļ���Ϣfeatures�棨features��ͶӰ�ں��������������ع�����ȷ���ܸ���һЩ��*/
			if (reConstructFeature(flatTerrainGd.get(), gg.get(), terrainBoundLoop.get()))
			{
				std::string name;
				if (g->getName().empty()==false)  name = g->getName(); else name = f->getName();
				gg->setName(name);
				
				gmList.push_back(gg);
			}
		}	//geode
		
		/* ����ÿ��geode������geometry, �����Χ���㼯�ϡ�*/
		for( std::vector<osg::ref_ptr<osg::Geometry> >::iterator qm= gmList.begin(); qm != gmList.end(); qm++ )
		{
			/* �����б߽�ĵ������������߽���ϸ������������ߵ��ཻ��������߽������Ч�߽�㡣*/
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

			/* �����߻��ļ��� */
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
				/* �Ȱѵ�ǰperlist�߻������ǻ���һ���棬Ϊ�˸��߽�μ������ཻ��� */				
				osg::ref_ptr<osg::Geometry> fg = ConvertFromLINELOOPToTranglesFace(*q);
				if(fg == NULL)	continue;
				if((qm->get()->getName().substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING))
				{
					lakeGd->addDrawable(fg);
				}

				/* �鿴 terrainBorderVertices �����б߽���Ƿ��� ��ǰperlist�߻��� �ཻ������ཻ����ô newIdx �� = 1������ = 0 */
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
	
						/* ����һ���жϽ�����ϵ��� */
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
				
				/* �鿴 terrainInnerVertices �������ڲ����Ƿ��� ��ǰperlist�߻��� �ཻ������ཻ����ô AllInnerVertexIntersection �� = true������ = false */
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
						
						/* ����ڲ�����perlist�߻����ཻ����ô����ڲ��������Ϊ������������*/
						if(isDetail == false)	innerVerticesInterWithCulture.push_back(qq->pVF);
					}
				}
				if(TotalInnerIntersectVerts == terrainInnerVertices.size()) AllInnerVertexIntersection = true;
				if(TotalInnerIntersectVerts == 0) AllInnerVertexNoIntersection = true;
	
				/* ���Ӽ�����Ϣ��ؾ��߽�Ľ���� */
				needDoSkirt = true;
				osg::ref_ptr<osg::Vec3Array> pmForCulture = new osg::Vec3Array;
				addIntersectVertex(flatTerrainGd, currPerlistGd, *q, pmForCulture, gmName);
	
				for(std::vector< osg::ref_ptr<osg::Vec3Array> >::iterator mq = multiDc.begin(); mq != multiDc.end(); mq++)
				{
					if (mq->get()->size()>=3) 
					{
						// dc ��������ϸ���������
						std::string name = qm->get()->getName();

						/* ��·, ���ȴ����� */
						if(name.substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)	roadList.push_back(*mq);
						if(name.substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)	hasLakes = true;

						osg::ref_ptr<ArealConstraint> dc = new ArealConstraint;
						if (isDetailedtxt(name))
						{
							hasDetailed = true;
						}
						else
						{
							/* �����Χ�����ϣ���ר��������ȹ���õġ�*/
							if(needDoSkirt)
							{
								/* ����ж�̬�ݣ�����Ҫ��ȹ�ߡ�*/
								if(name.substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) != DEFINE_DYNAMICGRASS_STRING)
								{
									perimeterList.push_back(*mq);
								}
							}
						}

						/* dc �����������Ļ���Ϣ���� */
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

	/* �����еĵ�����Ķ��㣬���պ������һ�¸߶� */
	TuneElevationOfLakeSurface(lakeGd, coords);

	/* ���Ʊ߽��ļ��� */
	/* ��ǰ��ѭ�����棬�Ѿ�������еġ��Ļ���Ϣ����α߽罻�㡱������allBorderVertices��*/
	/* ����һ���������е��α���ı߽�㣬Ҳ���뵽 allBorderVertices */
	for( std::vector<FindValidVertex>::iterator q= terrainBorderVertices.begin(); q != terrainBorderVertices.end(); q++ )
		allBorderVertices.push_back(*q);

	/* �Ȱ����еĵ�����Ķ��㣬�����е�totalVertices�������� */
	getTotalFaceVertex();

	/* ��·���浽ϸ�������� */
	if(roadList.size() > 0)
	{
		for(std::vector< osg::ref_ptr<osg::Vec3Array> >::iterator q = roadList.begin(); q != roadList.end(); q++)
		{
			/*   ��һ��    */
			/* �ѵ�ǰperlist�߻������ǻ���һ���� */				
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

			/* ������������ */
			osg::ref_ptr<osg::Vec3Array> validCoords = dynamic_cast<osg::Vec3Array *>(gm->getVertexArray());

			/* ��ƽ���ӳ��ؿռ�� */
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
				
			/*   �ڶ���    */
			/* �ѵ�ǰperlist�߻�������һ��ȹ�� */
			osg::ref_ptr<osg::Geometry> gmSkirt = ConvertFromLINELOOPToSkirtFace(*q);
			if(gmSkirt == NULL) continue;

			/* ������������ */
			osg::ref_ptr<osg::Vec3Array> skirtCoords = dynamic_cast<osg::Vec3Array *>(gmSkirt->getVertexArray());

			/* ��ƽ���ӳ��ؿռ�� */
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


	/* ��ϸ�������dc�������ó��� */
	MoveOutDetailTxtDc();

	/* ����߽�㣬�Ѵ��ⱻ�Ļ���Ϣ���ǵĵ���뵽��Ч�㼯��invalidVertices���� */
	int numInvalidBorderVertex = 0;
	bool allBorderInvalid = false;

#if 0	//���������������ȷ�Ļ�����γ���ʵ����Ҳ���Բ�Ҫ��
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

	/* ��ͬ���Ķ��dc�ϲ���һ����ϳ�һ��dc */
	/* �����о�����ͬ���Ƶ�dc�ϲ���Ϊһ���������㸲����� */
	MergeDclistStepOne(&dclist);

	/* ���dclist����ֻ������dc, ������dcһģһ������ô����������һ����*/
	CheckDclistAndMerge(&dclist);

	//���û����Ч��Լ�����ɣ������˳�
	if (dclist.empty() && (dcDetail.empty() && dcLakes.empty() ))
	{
		/* û����Ч������Լ�����ɵĻ���...... */
		if (IntersectCt < faceVertice.size())			// < 3
		{
			/* ��ǰ���������Ȼû�б��и�, �����Ƿ񱻵����˸߶�. */
			if(validCSM.size())
			{
				/* ��Ҫ�����޸ĺ��faceVertice����ԭ��coords�����Ӧ�Ķ���ֵ���������ֶϲ� */
				/* ��Ϊ���ж�����ţ�����ֱ�Ӱ������ָʾ�޸ļ��� */
				for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
				{
					int idx = f->idx;	terrainVertices[idx].pV = f->pV;
				}
				int index = 0;
				for( osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
				{
					*r = (terrainVertices[index++].pV);
				}

				/* ���˸߶ȵ�ҲҪ����д�ؾ� */
				CurrentTerrainIsCut = 1;
			}
			
#ifdef	MAKE_HIGH_SKIRT
			/* ����ȹ���б��ӳ� ����һ��ȹ�ߵĶ��㡣2013-01-04 */
			if(1)
			{
				/* ��Ϊ���ж�����ţ�����ֱ�Ӱ������ָʾ�޸ļ��� */
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

	/* ����ؾ��б仯������ȫ��ͨѶ������*/
	CurrentTerrainIsCut = 1;

	//�������1�����ϵĵ��feature���ཻ����Ҫ���ؾ�
	int VertexNoIntersect = abs((int)faceVertice.size() - IntersectCt );
	if(( VertexNoIntersect < 1) && (!hasDetailed) && (dclist.size() == 1))
		return NULL;

	/***********************************************************************/
	/********  ����ϸ��������  *********************************************/
	if((dcDetail.size() > 0) && (dclist.size() == 0))
	{
		/* ����������Ϣ���� */
		CreateFaceOfDetailTxtWhenNoCulture(&tcp);

		/* ����������һ��״�����Ļ���Ϣ��ϸ������ռ�����������Σ�������Ӧ����ǰ���˳���
		   ����Ҫ���ϸ������֮������˳�����������ȱ����һ��ϸ���������ˡ�
		 */ 
		if((VertexNoIntersect < 1 ) && (hasDetailed))
			return NULL;
		
		/* ��dcDetail��ֵ��dclist */
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcDetail.begin(); p != dcDetail.end(); p++ )
		{
			dclist.push_back(*p);
		} 
		dcDetail.clear();
	}

	/***********************************************************************/
	/********  ��������  *********************************************/
	if( (dcLakes.size() > 0) )
	{
		/* ���������������������˵�������ؿ鶼�Ǻ��棬Ҫ��ǰ�˳���*/
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
		
		/* ��Lake��Island��dc�ϲ���Ϊһ���������㸲����� */
		MergeLakeAndIsland(&dclist);

		/* ���ɺ��� */
		CreateFaceOfLake(&tcp);
		
		/* ����е��Ļ����ѵ���dclist����ȥ�� */
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

		/* ��dcLakes ��ֵ�� dclist */
		for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dcLakes.begin(); p != dcLakes.end(); p++ )
		{
			dclist.push_back(*p);
		} 
		dcLakes.clear();
	}


	/***********************************************************************/
	/********  ������̬�ݵ���  *********************************************/
	if((dclist.size() != 0) )
	{

		/* ���dclist�����Ƿ���DynamicGrass�����û�о��˳���*/
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
				/* ���ɶ�̬�ݵ��� */
				CreateFaceOfDynamicGrass(&tcp);
			}
	
			/* ����̬��dc��dclist��ȥ�� */
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


	//����һ���µĵ㼯�ϣ�ֻ�������ϵĵ�
	osg::ref_ptr<osg::Vec3Array> geVx = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator f= faceVertice.begin(); f != faceVertice.end(); f++ )
	{
#if 0			/* 2013-2-26 Galen : ���Է�����ǰ���Ѿ���������������ˡ�*/
		/* 
		   ������ϵĵ��Լ���ϵĵ����ܽ������ϵ������Ͳ��ӽ��������ǻ��������������ǻ���������ΪԼ���ϵĵ�����ϵĵ���������ʧԼ���ϵĵ㲢�����ն���������֣�
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

	//���Ĳ���1������ԭ��ȥ����ȹ�ߵĵ㣬��ȷ�����е����ȷ���

	/* �����ǻ�������е㼯�϶�����validVertices���� */
	/* ����validVertices�������е㣬��ͶӰ��ӳ�䵽�ռ�� */
	j = 0;
	validVertices.clear();	validElement.clear();
	for( osg::Vec3Array::iterator r = dcoords->begin(); r != dcoords->end(); r++ )
	{
		FindValidVertex vVert;
		osg::Vec3 pSpace = getSpaceVertex(*r);
		vVert.pV = pSpace;	vVert.pVF = *r; vVert.idx = j; vVert.newIdx = j; j++;
		validVertices.push_back(vVert);
	}

#if 0		//��������������ȷ�Ļ�����γ�����Բ�Ҫ��
	// �����ȥ���Ļ���Ϣ�ڲ�����Ч�棬��ôҪ�ѵ����ڲ����Ļ���Ϣ�ཻ�ĵ㶼��Ϊ��Ч�㡣
	for( std::vector<osg::Vec3>::iterator p= innerVerticesInterWithCulture.begin(); p != innerVerticesInterWithCulture.end(); p++ )
	{
		invalidVertices.push_back( *p );
	}
#endif

	//�ҳ����е���Ч�档����� validElement ����
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

	//����һ���µĵ����У�����������Ч��
	osg::ref_ptr<osg::Vec3Array> validCoords = new osg::Vec3Array;
	for( std::vector<FindValidVertex>::iterator pp= validVertices.begin(); pp != validVertices.end(); pp++ )
		validCoords->push_back(pp->pV);
	
	//��ֻ�����Ϸ������ε�validElement�����Ԫ�ؼӽ���
	osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(GL_TRIANGLES);
	for(std::vector<unsigned int>::iterator p = validElement.begin(); p!=validElement.end(); p++)
		dui->push_back(*p);

	if (dui->empty())	return NULL;

#ifdef ADD_SKIRT
	/* ����ȹ�� */
	CreateSkirt(validCoords, dui, 1);
#endif

	/* ������������и��ļ����塣*/
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


	/* �����ϸ�������ͼ�ϸ������ */
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

/* ����ϸ��������Ʋ� */
void ChangeCoordSystem(osg::Geode* gd, const osg::Matrix mt);
void AddDetailedTxtIntoTerrain(osg::MatrixTransform* mt, bool bGroup)
{
	/* ��·�����ó�����*/
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
	
	/* �����·���ȴ���·����·Ҳ���ɶ���������ƣ�������ʾԭ������ҹ����ʾ��������*/
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


	/* �Ѻ�Ҳ�����ó�����*/
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
	
	/* ����к����ȴ�������Ѻ���*/
#if 1
	/* IVE�ļ�������ʱ��������� */
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

	/* �ٴ���ϸ������ */
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
						
						//ֱ�ӷ���root����
						root->addChild(gd);

						//���н�㻻����������ϵ 
						ChangeCoordSystem(gd, mt->getMatrix());

						/* ���������и�*/
						process(root, false);

						/* ���н�㻹ԭΪ�ֲ�����ϵ */
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
				
				//ֱ�ӷ���root����
				root->addChild(gd);

				//���н�㻻����������ϵ 
				ChangeCoordSystem(gd, mt->getMatrix());

				/* ���������и�*/
				process(root, false);

				/* ���н�㻹ԭΪ�ֲ�����ϵ */
				ChangeCoordSystem(gd, osg::Matrix::inverse(mt->getMatrix()));
			}
			AddDetailedTxtIntoTerrain(mt, false);
			AddDynamicGrassIntoTerrain(mt, false);
		}
	}
	return 1;
}

/* ��������ϵ��ԭΪ��������ϵ����������ܰ���������Ʋ㣩*/
void RestoreToWorldCoordinate(osg::MatrixTransform* mt)
{
	if(mt)
	{
		for(unsigned int k=0; k<mt->getNumChildren();++k)
		{
			osg::Geode* gd = dynamic_cast<osg::Geode*>(mt->getChild(k));
			if(gd != NULL)
			{		
				/* ���н�㻹ԭΪ�ֲ�����ϵ */
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
							/* ���н�㻹ԭΪ�ֲ�����ϵ */
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
								/* ���н�㻹ԭΪ�ֲ�����ϵ */
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
						
						//���н�㻻����������ϵ 
						ChangeCoordSystem(gd, mt->getMatrix());
	
						//ֱ�ӷ���root����
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
						
						//���н�㻻����������ϵ 
						ChangeCoordSystem(gd, mt->getMatrix());

						//ֱ�ӷ���root����
						root->addChild(gd);
	
					}
				}
			}
		}

		/* ���������и�*/
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

/* ��ȡ��ǰIVE�ؾ��ļ���LOD���� */
void GetCurrentTerrainLodLevel(char *fnIve)
{
	/* �������ĵ����ļ��ǲ��Ƕ��������ļ������õ�����·���ļ��ļ���fn. */
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
	
	/* ���»��߽����Ƿֿ���*/
	splitMaterialName(fn);
	if(materialParts.size() > 2)
	{
		char levelStr[2]; levelStr[0] = materialParts[1][1]; levelStr[1] = 0;
		CurrentTerrainLodLevel = atoi(levelStr);
	}
	
}

int UpdateIVEByBTG(osg::ref_ptr<osg::Group> btg, char *fnIve)
{
	/* �и�ؾ�ǰ����λȫ��ͨѶ������*/
	CurrentTerrainIsCut = 0;
	UpdateTheElevationInTerrain = 0;
	allCultures = btg;

	/* ��ȡ��ǰIVE�ؾ��ļ���LOD���� */
	GetCurrentTerrainLodLevel(fnIve);

	/* ��ʼ��������Ϣ�б� */
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
