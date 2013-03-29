#ifndef OSGTDS_TERRAIN_PROCESSING_UTILS_DEF
#define OSGTDS_TERRAIN_PROCESSING_UTILS_DEF 1

#include <stack>

#include <osg/Vec3>
#include <osg/MatrixTransform>
#include <osg/Geode>

#include <osg/NodeVisitor>

struct VerIdx		//单个奇点
{
	osg::Vec3 Ver;
	int		  idx;
	int			times;
};

struct VerPair		//奇点对
{
	osg::Vec3 Ver;
	int		  start;
	int		  end;
};

/* 合法顶点数据结构 */
typedef struct _PROJECTIONMAP_
{
	osg::Vec3 pV;
	osg::Vec3 pVF;
	osg::Vec2 pT;
	double	  lon;
	double	  lat;
	double	  elev;
}ProjectionMap;

class TerrainProcessingUtils {

    public:
        // Test if a point is within a hull
        static bool isWithin( osg::Vec3 pt, osg::Vec3Array *hull );

        // Create a bottomHull from a list of vertices under Node
        static osg::ref_ptr<osg::Vec3Array> bottomHull( osg::Node &node );

        // Create a convex hull  from a list of vertices (Used by bottomHull)
        static osg::ref_ptr<osg::Vec3Array> convexHull( const osg::Vec3Array &v );

        // Compute a convex periphery (not perimeter).
        static osg::ref_ptr<osg::Vec3Array> computePeriphery( const osg::Vec3Array &coords);

        // Get the (tight) perimeter around the geometry g, returned in coords
        static void computePerimeter( osg::Geometry *g, osg::Vec3Array &coords);
		static void computePerimeterAll( osg::Geometry *g, std::vector<osg::ref_ptr<osg::Vec3Array>> &perimeterList );

        // Find the intersect point i, where point v intersects the terrain db
        static bool getIntersect( osg::Geode *db, const osg::Vec3 &p, osg::Vec3 &i); 

        // Find the intersect point i, where point v intersects the terrain db
        static bool getIntersect( osg::Group *db, const osg::Vec3 &p, osg::Vec3 &i); 

        // Reduce array vector in_v to a unique set of points.
        static osg::ref_ptr<osg::Vec3Array> uniqueCoords( osg::Vec3Array *in_v );

        // returns the Primitive type under node and fills in points with the coordinates
        static GLenum identifyPrimitiveAndGetPoints(osg::Node &node, osg::ref_ptr<osg::Vec3Array> &points);

        class IDPrimAndGetPoints: public osg::NodeVisitor
        {
            public:
                IDPrimAndGetPoints(osg::ref_ptr<osg::Vec3Array> &);
                virtual void apply(osg::Geode& geode);
                GLenum getPrimType();

            private:

                osg::ref_ptr<osg::Vec3Array> &_points;
                GLenum _primType;
        };

        static bool hasBottomHull( osg::Node &node );
        class HasBottomHull : public osg::NodeVisitor
        {
            public:
                HasBottomHull();
                virtual void apply(osg::Geode& geode);
                bool hasBottomHull();

            public:
                bool _flag;
                bool _flagSet;
        };


        class FindLowestZ : public osg::NodeVisitor
        {
            public:
                FindLowestZ();
                virtual void apply(osg::MatrixTransform& tx);
                virtual void apply(osg::Geode& geode);
                double getLowestZ();

            private:
                double _lowestZ;
                osg::Matrix _mat;
                std::stack<osg::Matrix> _matStack;
        };

        class GetBottomCoords : public osg::NodeVisitor
        {
            public:
                GetBottomCoords(double z, double epsilon=0.025);
                virtual void apply(osg::MatrixTransform& tx);
                virtual void apply(osg::Geode& geode);
                const osg::Vec3Array &getVertexArray();
        
            private:
                double _z;
                double _epsilon;
                osg::ref_ptr<osg::Vec3Array> _varray;
                osg::Matrix _mat;
                std::stack<osg::Matrix> _matStack;
        };


        class FindPeriphery : public osg::NodeVisitor
        {
            public:
                FindPeriphery();
                virtual void apply(osg::MatrixTransform& tx);
                virtual void apply(osg::Geode& geode);
                osg::ref_ptr<osg::Vec3Array> getPeriphery();

            private:
                double _z;
                osg::ref_ptr<osg::Vec3Array> _varray;
                osg::ref_ptr<osg::Vec3Array> _periph;
                osg::Matrix _mat;
                std::stack<osg::Matrix> _matStack;
        };

        struct Edge
        {
            osg::Vec3 A, B;
            bool keep;
        
            Edge(const Edge  &e );//: A(e.A), B(e.B), keep(true) {}
            Edge(const osg::Vec3 &a, const osg::Vec3 &b );//: A(a), B(b) {}
        
            bool operator == (Edge &rhs );
                /*
            {
                if( (A == rhs.A && B == rhs.B) || (A == rhs.B && B == rhs.A))
                    return true;
                else
                    return false;
            }
            */
        
            bool has( osg::Vec3 v );
            /*
            {
                if( v == A || v == B ) return true;
                return false;
            }
            */

			double length();
        };

        struct MakeEdgeList
        {
            public:
                MakeEdgeList();/*:
                    edges(0L) {}*/
        
                void setEdgeVector( std::vector<Edge> *e );
                /*
                {
                    edges = e;
                }
                */
        
                void operator() (const osg::Vec3& v1,const osg::Vec3& v2,const osg::Vec3& v3, bool flag) const;
                /*
                {
                    if( edges != 0L )
                    {
                        edges->push_back(Edge( v1, v2 ));
                        edges->push_back(Edge( v2, v3 ));
                        edges->push_back(Edge( v3, v1 ));
                    }
                }
                */
            private:
        
                std::vector<Edge> * edges;
        };

        // 判断当前线环是否合法
		static Edge* CheckLINELOOPValibility(osg::ref_ptr<osg::Vec3Array> coords);
};



#endif

