#ifndef OSGTDS_CORRECTION_SPACE_MESH_H
#define OSGTDS_CORRECTION_SPACE_MESH_H

#include <osg/Referenced>
#include <osg/Array>
#include <osg/Vec3>

typedef struct _LIMITVALUE_
{
	double maxLongi;
	double maxLati;
	double minLongi;
	double minLati;
}LimitValue;

class CorrectionSpaceMesh : public osg::Referenced
{
    public:
        CorrectionSpaceMesh( osg::Group *targets, double K);

        void correctPoint( osg::Vec3 &p, bool print=false );
        osg::Vec3 getZCorrection( const osg::Vec3 &v );

        void  calcCSMGeode();
        osg::ref_ptr<osg::Geode> getCSMGeode() {	return csmGd; }

        void setK(double K)  { _K = K; }
        double getK() { return _K; }

        const std::vector< osg::ref_ptr<osg::Vec3Array> > &getBottomHulls() const { return  _bottomHulls; }
        const osg::ref_ptr<osg::Vec3Array>& getNonPolyTargetVerts() const { return  _nonPolyTargetVerts; }

				/* 包围球的球心和半径 */
        osg::Vec3 center;
        float radius;
        osg::ref_ptr<osg::Geode> csmGd;

    private:
        bool _valid;
        double _K;

        class Triangle {
            public:
                Triangle( osg::Vec3 _a, osg::Vec3 _b, osg::Vec3 _c );

                void setK( double K ) { _K = K; }
                double getK() { return _K; }

                void getABC( osg::Vec3 &A, osg::Vec3 &B, osg::Vec3 &C );
        
                void print();
        
                bool within( double m, double n );
                bool within( const osg::Vec3d &v );

                // Height Correction Function
                // C(x,y) = w1*f1(x,y) + w2*f2(x,y) + w3*f3(x,y)
                double C( const osg::Vec3d &p );
        
            private:
                osg::Vec3 a, b, c;			// a,b,c为三角形三顶点的直角坐标系的值
            	osg::Vec3d x, y, z;			// x,y,z为三角形三顶点的经纬高坐标系的值，x=a; y=b; z=c;
                
                double _K;
        
                double _length2D( osg::Vec3d A, osg::Vec3d B );
                bool   _eqline( osg::Vec3d A, osg::Vec3d B, double &m, double &d );
                double _det( osg::Vec3d A, osg::Vec3d B, osg::Vec3d P );
                void   _intersect( osg::Vec3d A, osg::Vec3d B,
                                 osg::Vec3d C, osg::Vec3d D,
                                 osg::Vec3d &I );
        
                double f(osg::Vec3d p, double m, double n );

        };

        std::vector<Triangle> _triangles;
        std::vector< osg::ref_ptr<osg::Vec3Array> > _bottomHulls; 
        osg::ref_ptr<osg::Vec3Array> _nonPolyTargetVerts;

};
#endif
