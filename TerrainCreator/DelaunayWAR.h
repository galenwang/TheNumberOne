#ifndef OSGTDS_DELAUNAY_WAR_DEF
#define OSGTDS_DELAUNAY_WAR_DEF

#include <set>
#include <vector>
#include <osgUtil/DelaunayTriangulator>

class DelaunayWAR: public osg::Referenced
{
    public:
        DelaunayWAR( osg::Vec3Array * );
        bool triangulate();
        const osg::DrawElementsUInt *getTriangles() const;
        osg::DrawElementsUInt *getTriangles();

    private:

        osg::ref_ptr<osg::Vec3Array> _coords;
        osg::ref_ptr<osgUtil::DelaunayTriangulator> _dt;

        struct PairPts 
        {
            osg::Vec3 A, B;
            PairPts( osg::Vec3 a, osg::Vec3 b ): A(a), B(b) {}
        };

        osg::ref_ptr<osg::Vec3Array> _periphery( const osg::Vec3Array &v );

        //bool _apply( osg::Vec3Array &coords, osgUtil::DelaunayTriangulator &dt );
        bool _apply();
        
        
        bool _pairInTriangles( 
                const osg::Vec3 &A, const osg::Vec3 &B,  
                const osg::Vec3Array &coords, 
                osg::DrawElementsUInt *t );
        
        /*
        std::set<osg::Vec3> _vec3ListOfTrianglesPtBelongsTo( osg::Vec3 P, 
                const osg::Vec3Array &coords, 
                osg::DrawElementsUInt *t );
                */
};

#endif
