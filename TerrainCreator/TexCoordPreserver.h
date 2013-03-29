#ifndef OSGTDS_TEXCOORD_PRESERVER_DEF
#define OSGTDS_TEXCOORD_PRESERVER_DEF

#include <osg/Geometry>
#include <osg/Vec3>
#include <osg/Vec2>
#include <osgUtil/DelaunayTriangulator>

#define MINDISTTOMERGEVERTEX 1.0

/* 合法顶点数据结构 */
struct FindValidVertex
{
	osg::Vec3 pV;
	osg::Vec3 pVF;
	osg::Vec3 pN;
	osg::Vec2 pT;
	int idx;
	int newIdx;
	int isValid;

	/* 赋值 */
	void set(FindValidVertex &rhs )
	{
		pV      = rhs.pV;
		pVF     = rhs.pVF;
		pN      = rhs.pN;
		pT      = rhs.pT;
		idx     = rhs.idx;
		newIdx  = rhs.newIdx;
		isValid = rhs.isValid;
	};
};

typedef struct _FINDBORDERVERTER_
{
	FindValidVertex fvv;
	int onDetailed;
	int onCulture;
}FindBorderVertex;

namespace osgTDS {

class TexCoordPreserver
{
	public:
		TexCoordPreserver(std::vector<FindValidVertex> *g);
		osg::ref_ptr<osg::Vec2Array> getTexCoords( const osg::Vec3Array &coords );
		osg::ref_ptr<osg::Vec3Array> getNormals( const osg::Vec3Array &coords );

	private:
		std::vector<FindValidVertex> *preserver;
		FindValidVertex Vorigin, Vxmax, Vymax;
		double xLen, yLen;
};

}



/* 克隆德罗尼三角化类 */
class TerrainTriangulator: public osg::Referenced 
{
public:
    TerrainTriangulator();
    explicit TerrainTriangulator(osg::Vec3Array *points, osg::Vec3Array *normals = 0);
    TerrainTriangulator(const TerrainTriangulator &copy, const osg::CopyOp &copyop = osg::CopyOp::SHALLOW_COPY);

    typedef std::vector< osg::ref_ptr<osgUtil::DelaunayConstraint> > linelist;

    /** Set the input point array. */
    inline void setInputPointArray(osg::Vec3Array* points) { points_ = points; }

    /** Get the const input point array. */
    inline const osg::Vec3Array* getInputPointArray() const {  return points_.get(); }

    /** Get the input point array. */
    inline osg::Vec3Array* getInputPointArray() {  return points_.get(); }


    /** Set the output normal array (optional). */
    inline void setOutputNormalArray(osg::Vec3Array* normals) { normals_ = normals; }

    /** Get the const output normal array (optional). */
    inline const osg::Vec3Array *getOutputNormalArray() const { return normals_.get(); }

    /** Get the output normal array (optional). */
    inline osg::Vec3Array *getOutputNormalArray() { return normals_.get(); }


    /** Add an input constraint loop.
     ** the edges of the loop will constrain the triangulation.
     ** if remove!=0, the internal triangles of the constraint will be removed;
     ** the user may the replace the constraint line with an equivalent geometry.
     ** GWM July 2005 */
    void addInputConstraint(osgUtil::DelaunayConstraint *dc) { constraint_lines.push_back(dc); }


    /** Start triangulation. */
    bool triangulate();

    /** Get the generated primitive (call triangulate() first). */
    inline const osg::DrawElementsUInt *getTriangles() const { return prim_tris_.get(); }

    /** Get the generated primitive (call triangulate() first). */
    inline osg::DrawElementsUInt *getTriangles() { return prim_tris_.get(); }
    
    /** remove the triangles internal to the constraint loops.
     * (Line strips cannot remove any internal triangles). */
    void removeInternalTriangles(osgUtil::DelaunayConstraint *constraint);


protected:
    virtual ~TerrainTriangulator();
    TerrainTriangulator &operator=(const TerrainTriangulator &) { return *this; }
    int getindex(const osg::Vec3 &pt,const osg::Vec3Array *points);

private:
    osg::ref_ptr<osg::Vec3Array> points_;
    osg::ref_ptr<osg::Vec3Array> normals_;
    osg::ref_ptr<osg::DrawElementsUInt> prim_tris_;

    // GWM these lines provide required edges in the triangulated shape.
    linelist constraint_lines;

    void _uniqueifyPoints();
};

#endif
