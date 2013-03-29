#include <osg/Geometry>
#include <osg/LineSegment>
#include <osgUtil/IntersectVisitor>
#include "TerrainProcessingUtils.h"
#include <algorithm>
#include <osg/TriangleFunctor>

#if defined( WIN32) && !defined(__CYGWIN__) && !defined( M_PI )
    // Grrr.
#   define M_PI 3.1415926535897932384626433832795
#endif


osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::uniqueCoords( osg::Vec3Array *in_v )
{
    std::sort( in_v->begin(), in_v->end() );

    osg::ref_ptr<osg::Vec3Array> nv = new osg::Vec3Array;
    osg::Vec3 v;
    for( osg::Vec3Array::iterator p = in_v->begin(); p != in_v->end(); p++ )
    {
        if( p == in_v->begin() )
        {
            v = *p;
            nv->push_back(v);
            continue;
        }

        if( *p == v )
            continue;

        v = *p;
        nv->push_back( v );
    }
    return nv;
}

bool TerrainProcessingUtils::isWithin( osg::Vec3 pt, osg::Vec3Array *hull )
{
    if( hull->size() == 1 )
    {
        if( pt == hull->front() )
            return true;
        else
            return false;
    }

    std::vector<osg::Vec3> _hull;
    for( osg::Vec3Array::iterator p = hull->begin(); p != hull->end(); p++ )
        _hull.push_back( *p);

    osg::Vec3 A, B;
    unsigned int n = _hull.size();
    for( unsigned int i = 0; i <= n; i++ )
    {
        if( i == 0 )
        {
            A = _hull[i];
            continue;
        }
        B = _hull[i%n];

        osg::Vec3 V0 = (B - A);
        osg::Vec3 V1 = (pt - A); 
        float d = (V0[0] * V1[1]) - (V0[1] * V1[0]);
        if( d < 0.0 )
            return false;

        A = B;
    }
    return true;
}

osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::bottomHull( osg::Node &node )
{
    // Find the coordinate in the model with the lowest altitude (Z)
    FindLowestZ flz;
    node.accept( flz );
    double Z = flz.getLowestZ();

    // Get all coordinates within an epsilon of Z
    GetBottomCoords gbc( Z );
    node.accept( gbc );

    osg::ref_ptr<osg::Vec3Array> newPoints; 

#if 0
    if( gbc.getVertexArray().size() == 1 ) // we have a point
    {
        //newPoints = new osg::Vec3Array;
        //newPoints->push_back( gbc.getVertexArray().front() );
        hull.push_back( gbc.getVertexArray().front() );
    }
    else
#endif
    newPoints  = convexHull( gbc.getVertexArray() );

    return newPoints;
}

bool TerrainProcessingUtils::getIntersect( osg::Geode *terrain, const osg::Vec3 &p, osg::Vec3 &i )
{
    // 10000.0 is arbitrary
    osg::Vec3 p0 = p + osg::Vec3(0,0,10000.0);
    osg::Vec3 p1 = p + osg::Vec3(0,0,-10000.0);

    osg::ref_ptr<osg::LineSegment> seg = new osg::LineSegment(p0,p1);

    osgUtil::IntersectVisitor iv;
    iv.addLineSegment( seg.get() );
    if (!seg->valid())
    {
        int i = 0;
	    return false;
    }
    else
    {
	    terrain->accept( iv );
	
	    if( iv.hits() )
	    {
	        osgUtil::IntersectVisitor::HitList& hitList = iv.getHitList(seg.get());
	        for(osgUtil::IntersectVisitor::HitList::iterator hitr=hitList.begin();
	            hitr!=hitList.end(); ++hitr)
	        {
	            i = hitr->getWorldIntersectPoint();
	            return true;
	        }
	    }
	}
    return false;
}

bool TerrainProcessingUtils::getIntersect( osg::Group *terrain, const osg::Vec3 &p, osg::Vec3 &i )
{
    // 10000.0 is arbitrary
    osg::Vec3 p0 = p + osg::Vec3(0,0,10000.0);
    osg::Vec3 p1 = p + osg::Vec3(0,0,-10000.0);

    osg::ref_ptr<osg::LineSegment> seg = new osg::LineSegment(p0,p1);

    osgUtil::IntersectVisitor iv;
    iv.addLineSegment( seg.get() );
    if (!seg->valid())
    {
        int i = 0;
	    return false;
    }
    else
    {
	    terrain->accept( iv );
	
	    if( iv.hits() )
	    {
	        osgUtil::IntersectVisitor::HitList& hitList = iv.getHitList(seg.get());
	        for(osgUtil::IntersectVisitor::HitList::iterator hitr=hitList.begin();
	            hitr!=hitList.end(); ++hitr)
	        {
	            i = hitr->getWorldIntersectPoint();
	            return true;
	        }
	    }
	}
    return false;
}

TerrainProcessingUtils::FindLowestZ::FindLowestZ() : osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN )
{
    _lowestZ = DBL_MAX;
    _mat = osg::Matrix::identity();
}

void TerrainProcessingUtils::FindLowestZ::apply(osg::MatrixTransform& tx)           
{ 
    _matStack.push( _mat );

    _mat = _mat * tx.getMatrix();

    traverse( tx );

    _mat = _matStack.top();
    _matStack.pop();
}


void TerrainProcessingUtils::FindLowestZ::apply(osg::Geode& geode)                     
{ 
    for( unsigned int i = 0; i < geode.getNumDrawables(); i++ )
    {
        osg::Geometry *geom = dynamic_cast<osg::Geometry *>(geode.getDrawable(i));

        if( geom )
        {
            osg::Vec3Array *varray = dynamic_cast<osg::Vec3Array *>(geom->getVertexArray());

            if( varray != 0L )
            {
                for( unsigned int j = 0; j < varray->size(); j++ )
                {
                    osg::Vec3Array::iterator p;
                    for( p = varray->begin(); p != varray->end(); p++ )
                    {
                        osg::Vec3 v = osg::Vec3(*p) * _mat;
                        if( v[2] < _lowestZ )
                            _lowestZ = v[2];
                    }
                }
            }
        }
    } 
}

double TerrainProcessingUtils::FindLowestZ::getLowestZ() 
{ 
    return _lowestZ; 
}

TerrainProcessingUtils::GetBottomCoords::GetBottomCoords(double z, double epsilon) : 
    osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN ), 
    _z(z), 
    _epsilon(epsilon) 
{
    _varray = new osg::Vec3Array;
    _mat = osg::Matrix::identity();
}

void TerrainProcessingUtils::GetBottomCoords::apply(osg::MatrixTransform& tx)           
{ 
    _matStack.push( _mat );
    _mat = _mat * tx.getMatrix();
    traverse( tx );
    _mat = _matStack.top();
    _matStack.pop();
}

void TerrainProcessingUtils::GetBottomCoords::apply(osg::Geode& geode)                     
{ 
    for( unsigned int i = 0; i < geode.getNumDrawables(); i++ )
    {
        osg::Geometry *geom = dynamic_cast<osg::Geometry *>(geode.getDrawable(i));

        if( geom )
        {
            osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(geom->getVertexArray());

            if( coords != 0L )
            {
                osg::Vec3Array::iterator p;
                for( p= coords->begin(); p != coords->end(); p++)
                {
                    osg::Vec3 v =  osg::Vec3(*p) * _mat;
                    if( fabs(v[2] - _z) < _epsilon )
                    {
                        bool found = false;
                        for( osg::Vec3Array::iterator q = _varray->begin(); q != _varray->end(); q++ )
                        {
                            osg::Vec3 w(*q);
                            if( w == v )
                            {
                                found = true;
                                break;
                            }
                        }
                        if( !found )
                        {
                            _varray->push_back( v );
                        }
                    }
                }
            }
        }
    } 
}

const osg::Vec3Array &TerrainProcessingUtils::GetBottomCoords::getVertexArray() 
{ 
    return *(_varray.get()); 
}


TerrainProcessingUtils::FindPeriphery::FindPeriphery() : 
    osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN )
{
    _varray = new osg::Vec3Array;
    _mat = osg::Matrix::identity();
}

void TerrainProcessingUtils::FindPeriphery::apply(osg::MatrixTransform& tx)           
{ 
    _matStack.push( _mat );
    _mat = _mat * tx.getMatrix();
    traverse( tx );
    _mat = _matStack.top();
    _matStack.pop();
}

void TerrainProcessingUtils::FindPeriphery::apply(osg::Geode& geode)                     
{ 
    for( unsigned int i = 0; i < geode.getNumDrawables(); i++ )
    {
        osg::Geometry *geom = dynamic_cast<osg::Geometry *>(geode.getDrawable(i));

        if( geom )
        {
            osg::Vec3Array *varray = dynamic_cast<osg::Vec3Array *>(geom->getVertexArray());

            if( varray != 0L )
            {
                osg::Vec3Array::iterator p;
                for( p= varray->begin(); p != varray->end(); p++)
                {
                    osg::Vec3 v =  osg::Vec3(*p) * _mat;
                    _varray->push_back( v );
                }
            }
        }
    } 
}

osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::FindPeriphery::getPeriphery() 
{ 
    _periph = TerrainProcessingUtils::computePeriphery( *(_varray.get()) );
    return _periph; 
}

osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::computePeriphery( const osg::Vec3Array &v )
{
    // Allocate the new array for the perimeter points
    osg::ref_ptr<osg::Vec3Array> newPoints = new osg::Vec3Array;

    // Find the "center" 
    osg::Vec3Array::const_iterator p = v.begin();
    float minX = (*p)[0],
          maxX = (*p)[0],
          minY = (*p)[1],
          maxY = (*p)[1];
    p++;
    for( ; p != v.end(); p++ )
    {
        if( (*p)[0] < minX ) minX = (*p)[0];
        if( (*p)[0] > maxX ) maxX = (*p)[0];
        if( (*p)[1] < minY ) minY = (*p)[1];
        if( (*p)[1] > maxY ) maxY = (*p)[1];
     }
    
    osg::Vec3 C( (minX + maxX)/2, (minY + maxY)/2, 0 );

    //  Find the point furthest (F) from the center
    osg::Vec3 F;
    double dist = 0.0;
    for( p = v.begin(); p != v.end(); p++ )
    {
        osg::Vec3 pp((*p)[0], (*p)[1], 0 );
        double l = (pp - C).length();
        if( l > dist )
        {
            F = *p;
            dist = l;
        }
    }

    // A = first point
    // L = previous (last) point
    // N = new point
    osg::Vec3 A = F;
    osg::Vec3 L = C;
    osg::Vec3 N;

    // The furthest point from the center will be on the perimeter
    newPoints->push_back( A );

    while( true )
    {
        // a = This point less Last point defines an orientation around 'z'
        // aa = the angle of that orientation
        osg::Vec3 a = A - L;
        double aa = atan2f( a[1], a[0] );
        if( aa < 0.0 )
            aa += M_PI*2.0;

        // Find the vector formed by A and each point in the array that
        // has the smallest angle difference between aa and bb
        double smallestAngle = 2*M_PI;
        for( p = v.begin(); p != v.end(); p++ )
        {
            // B = each point in the array (except A)
            osg::Vec3 B = *p;
            if( B == A )
                continue;

            // b = the vector between AB
            // bb = defines the orientation around 'z' of AB
            osg::Vec3 b = B - A;
            double bb = atan2f( b[1], b[0] ) - aa;
            if( bb < 0.0 ) bb += M_PI*2.0;

            // Record the smallest angle and the point
            if( bb < smallestAngle )
            {
                smallestAngle = bb;
                N = B;
            }
        }
        // If N == F we've circumvented the perimeter
        if( N == F )
            break;

        L = A;
        A = N;

        newPoints->push_back( N );
    }


    // Now filter out the points that lie in the interim between
    // peripheral points.
    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;
    {
        int i = 0;
        osg::Vec3  A, B, L;
        double angle = 0.0;
        for( osg::Vec3Array::iterator p = newPoints->begin(); 
                p != newPoints->end(); p++ )
        {
            ++i;
            if( i == 1 )
            {
                A = *p;
                points->push_back( A );
                continue;
            }
            else if( i == 2 )
            {
                B = *p;
                osg::Vec3 a = B - A;
                angle = atan2f( a[1], a[0] );
                if( angle < 0.0 )
                    angle += M_PI*2.0;
                continue;
            }
            else
            {
                L = B;
                B = *p;
                osg::Vec3 a = B - A;
                double aa = atan2f( a[1], a[0] );
                if( aa < 0.0 )
                    aa += M_PI*2.0;
                if( aa > angle )
                {
                    points->push_back( L );
                    A = L;
                    osg::Vec3 a = B - A;
                    angle = atan2f( a[1], a[0] );
                    if( angle < 0.0 )
                        angle += M_PI*2.0;
                }
                
            }
        }
        // Don't forget the last point...
        if( i > 2 )
            points->push_back( B );
    }

    return points.get();
}

#if 0
osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::FindPeriphery::_periphery( const osg::Vec3Array &v )
{
    // Allocate the new array for the perimeter points
    osg::ref_ptr<osg::Vec3Array> newPoints = new osg::Vec3Array;

    // Find the "center" (average)
    osg::Vec3Array::const_iterator p = v.begin();
    float minX = (*p)[0],
          maxX = (*p)[0],
          minY = (*p)[1],
          maxY = (*p)[1];
    p++;
    for( ; p != v.end(); p++ )
    {
        if( (*p)[0] < minX ) minX = (*p)[0];
        if( (*p)[0] > maxX ) maxX = (*p)[0];
        if( (*p)[1] < minY ) minY = (*p)[1];
        if( (*p)[1] > maxY ) maxY = (*p)[1];
     }
    osg::Vec3 C( (minX + maxX)/2, (minY + maxY)/2, 0 );

    //  Find the point furthest (F) from the center
    osg::Vec3 F;
    double dist = 0.0;
    for( p = v.begin(); p != v.end(); p++ )
    {
        double l = (*p - C).length();
        if( l > dist )
        {
            F = *p;
            dist = l;
        }
    }

    // A = first point
    // L = last point
    // N = new point
    osg::Vec3 A = F;
    osg::Vec3 L = C;
    osg::Vec3 N;

    // The furthest point from the center will be on the perimeter
    newPoints->push_back( A );

    while( true )
    {
        // a = This point less Last point defines an orientation around 'z'
        // aa = the angle of that orientation
        osg::Vec3 a = A - L;
        double aa = atan2f( a[1], a[0] );
        if( aa < 0.0 )
            aa += M_PI*2.0;

        // Find the vector formed by A and each point in the array that
        // has the smallest angle difference between aa and bb
        double smallestAngle = 2*M_PI;
        for( p = v.begin(); p != v.end(); p++ )
        {
            // B = each point in the array (except A)
            osg::Vec3 B = *p;
            if( B == A )
                continue;

            // b = the vector between AB
            // bb = defines the orientation around 'z' of AB
            osg::Vec3 b = B - A;
            double bb = atan2f( b[1], b[0] ) - aa;
            if( bb < 0.0 ) bb += M_PI*2.0;

            // Record the smallest angle and the point
            if( bb < smallestAngle )
            {
                smallestAngle = bb;
                N = B;
            }
        }
        // If N == F we've circumvented the perimeter
        if( N == F )
            break;

        L = A;
        A = N;

        newPoints->push_back( N );
    }


    osg::ref_ptr<osg::Vec3Array> points = new osg::Vec3Array;
    {
        int i = 0;
        osg::Vec3  A, B, L;
        double angle = 0.0;
        for( osg::Vec3Array::iterator p = newPoints->begin(); 
                p != newPoints->end(); p++ )
        {
            ++i;
            if( i == 1 )
            {
                A = *p;
                points->push_back( A );
                continue;
            }
            else if( i == 2 )
            {
                B = *p;
                osg::Vec3 a = B - A;
                angle = atan2f( a[1], a[0] );
                if( angle < 0.0 )
                    angle += M_PI*2.0;
                continue;
            }
            else
            {
                L = B;
                B = *p;
                osg::Vec3 a = B - A;
                double aa = atan2f( a[1], a[0] );
                if( aa < 0.0 )
                    aa += M_PI*2.0;
                if( aa > angle )
                {
                    points->push_back( L );
                    A = L;
                    osg::Vec3 a = B - A;
                    angle = atan2f( a[1], a[0] );
                    if( angle < 0.0 )
                        angle += M_PI*2.0;
                }
                
            }
    
        }
    }

    return points.get();
}
#endif

osg::ref_ptr<osg::Vec3Array> TerrainProcessingUtils::convexHull( const osg::Vec3Array &v )
{
    // Allocate the new array for the perimeter points
    osg::ref_ptr<osg::Vec3Array> newPoints = new osg::Vec3Array;

    if( v.size() == 1 )
    {
        newPoints->push_back( v.front() );
        return newPoints.get();
    }

    // Find the "center" (average)
    osg::Vec3Array::const_iterator p = v.begin();
    float minX = (*p)[0], 
          maxX = (*p)[0], 
          minY = (*p)[1], 
          maxY = (*p)[1];
    p++;
    for( ; p != v.end(); p++ )
    {
        if( (*p)[0] < minX ) minX = (*p)[0];
        if( (*p)[0] > maxX ) maxX = (*p)[0];
        if( (*p)[1] < minY ) minY = (*p)[1];
        if( (*p)[1] > maxY ) maxY = (*p)[1];
     }
    osg::Vec3 C( (minX + maxX)/2, (minY + maxY)/2, 0 );

    //  Find the point furthest (F) from the center
    osg::Vec3 F;
    double dist = 0.0;
    for( p = v.begin(); p != v.end(); p++ )
    {
        double l = (*p - C).length();
        if( l > dist )
        {
            F = *p;
            dist = l;
        }
    }

    // A = first point
    // L = last point
    // N = new point
    osg::Vec3 A = F;
    osg::Vec3 L = C;
    osg::Vec3 N;

    // The furthest point from the center will be on the perimeter
    newPoints->push_back( A );

    while( true )
    {
        // a = This point less Last point defines an orientation around 'z'
        // aa = the angle of that orientation
        osg::Vec3 a = A - L;
        double aa = atan2f( a[1], a[0] );
        if( aa < 0.0 ) 
            aa += M_PI*2.0;

        // Find the vector formed by A and each point in the array that 
        // has the smallest angle difference between aa and bb
        double smallestAngle = 2*M_PI;
        for( p = v.begin(); p != v.end(); p++ )
        {
            // B = each point in the array (except A)
            osg::Vec3 B = *p;
            if( B == A )
                continue;

            // b = the vector between AB
            // bb = defines the orientation around 'z' of AB
            osg::Vec3 b = B - A;
            double bb = atan2f( b[1], b[0] ) - aa;
            if( bb < 0.0 ) bb += M_PI*2.0;

            // Record the smallest angle and the point
            if( bb  < smallestAngle )
            {
                smallestAngle = bb;
                N = B;
            }
        }

        // If N == F we've circumvented the perimeter
        if( N == F )
            break;

        L = A;
        A = N;

        newPoints->push_back( N );
    }

    return newPoints.get();
}

bool TerrainProcessingUtils::hasBottomHull( osg::Node &node )
{
    HasBottomHull hbh;
    node.accept(hbh);
    return hbh.hasBottomHull();
}

TerrainProcessingUtils::HasBottomHull::HasBottomHull() : 
    osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN ),
    _flag(false),
    _flagSet(false)
{}


bool TerrainProcessingUtils::HasBottomHull::hasBottomHull()
{
    return _flag;
}

void TerrainProcessingUtils::HasBottomHull::apply(osg::Geode& geode)
{
    for( unsigned int i = 0; i < geode.getNumDrawables(); i++ )
    {
        osg::Geometry *g = dynamic_cast<osg::Geometry *>(geode.getDrawable(i));
        if( g )
        {
            osg::PrimitiveSet  *pset = g->getPrimitiveSet(0);

            switch( pset->getMode() )
            {
                case osg::PrimitiveSet::POINTS:
                case osg::PrimitiveSet::LINES:
                case osg::PrimitiveSet::LINE_STRIP:
                case osg::PrimitiveSet::LINE_LOOP:
                    _flag = false;
                    _flagSet = true;
                    break;
                    
                case osg::PrimitiveSet::TRIANGLES:
                case osg::PrimitiveSet::TRIANGLE_STRIP:
                case osg::PrimitiveSet::TRIANGLE_FAN:
                case osg::PrimitiveSet::QUADS:
                case osg::PrimitiveSet::QUAD_STRIP: 
                case osg::PrimitiveSet::POLYGON:
                    if( !_flagSet )
                    {
                        _flag = true;
                        _flagSet = true;
                    }
                    break;
                default: break;
            }
        }
    }
}

GLenum TerrainProcessingUtils::identifyPrimitiveAndGetPoints(osg::Node &node, osg::ref_ptr<osg::Vec3Array> &points)
{
    IDPrimAndGetPoints idpagp(points);
    node.accept( idpagp );
    return idpagp.getPrimType();
}

void TerrainProcessingUtils::IDPrimAndGetPoints::apply(osg::Geode& geode)
{
    // Assume only one drawable
    osg::Geometry *g = dynamic_cast<osg::Geometry *>(geode.getDrawable(0));
    if( g )
    {
        osg::PrimitiveSet  *pset = g->getPrimitiveSet(0);
        _primType = pset->getMode(); 
        osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
        for( osg::Vec3Array::iterator p = coords->begin(); p != coords->end(); p++ )
            _points->push_back( *p );
    }
    // Don't traverse.  We're done.
}

TerrainProcessingUtils::IDPrimAndGetPoints::IDPrimAndGetPoints( osg::ref_ptr<osg::Vec3Array> &points):
    osg::NodeVisitor( osg::NodeVisitor::TRAVERSE_ALL_CHILDREN ),
    _points(points),
    _primType(0xFFFF)
{
}

GLenum TerrainProcessingUtils::IDPrimAndGetPoints::getPrimType()
{
    return _primType;
}

void TerrainProcessingUtils::computePerimeter( osg::Geometry *g, osg::Vec3Array &coords )
{
    std::vector<Edge>   edges;
    osg::TriangleFunctor<MakeEdgeList> mel;
    mel.setEdgeVector(&edges);
    g->accept( mel );

    for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
    {
        std::vector<Edge>::iterator q = p;
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

    std::vector<Edge> kept;
    for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
    {
        if( p->keep )
            kept.push_back( *p );
    }

	Edge e = kept.front();
	osg::Vec3 A = e.A;
	osg::Vec3 B = e.B;
	do {
		for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
		{
			if( p->has(B) && !(e == *p) )
			{
				coords.push_back( A );
				A = B;
				if( B == p->A )
					B = p->B;
				else
					B = p->A;

				e = *p;
				break;
			}
		}
	}while( !(e == kept.front()) );
	coords.push_back(A);
}


/* ���һ���߻���������һ���ߵĳ��ȣ����ڻ�ӽ�������������ߵĳ���֮�ͣ���ô����߻���Ч����ı�Ӧ��ɾ����*/
osg::Vec3 Zero(0.0, 0.0, 0.0);
TerrainProcessingUtils::Edge ret(Zero, Zero);
TerrainProcessingUtils::Edge* TerrainProcessingUtils::CheckLINELOOPValibility(osg::ref_ptr<osg::Vec3Array> coords)
{
	/* ��һ�����ҳ�����߻������б� */
	std::vector<Edge>   edges;
	for( osg::Vec3Array::iterator p = coords->begin(); p != coords->end(); p++ )
	{
		if( p == coords->end() - 1 ) break;
		osg::Vec3Array::iterator q = p + 1;
		Edge e(*p, *q);	edges.push_back(e);
	}
	
	/* �ڶ�������ÿ���ߵĳ��ȣ����ҳ����ĳ��ȡ�*/
	double maxLen = 0.0, otherLens = 0.0;
	for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		double currLen = p->length();
		if(currLen > maxLen)	{	maxLen = currLen;	ret = *p;	}
		otherLens += currLen;
	}
	
	/* �������ж��ܳ����ǲ�����󳤶ȵ����� */
	if((otherLens - maxLen * 2.0) < 5.0)	return &ret;
	else return NULL;
}

/* �ҳ����е���� */
void FindOutAllOddPoints( osg::Geometry *g, std::vector<TerrainProcessingUtils::Edge> *edges, std::vector<VerIdx> *oddpoints, int condition)
{
	oddpoints->clear();
	osg::ref_ptr<osg::Vec3Array> vx = dynamic_cast<osg::Vec3Array *> (g->getVertexArray());
	int index = 0;
	for( osg::Vec3Array::iterator r = vx->begin(); r != vx->end(); r++ )
	{
		int times = 0;
		for( std::vector<TerrainProcessingUtils::Edge>::iterator p = edges->begin(); p != edges->end(); p++ )
		{
    		if(*r == p->A)	times++;
    		if(*r == p->B)	times++;
		}
		
		/* ���ݲ�ͬ�����������ͬ���ж�������*/
		bool ifCond;
		if(condition == 0)	 ifCond = (times == 1);
		if(condition == 1)	 ifCond = (times >  2);
			
		if(ifCond)
		{
			int iFlag = 0;
			for(unsigned int j=0; j< oddpoints->size(); j++)
			{
				VerIdx curr = oddpoints->at(j);
				if(*r == curr.Ver) { iFlag = 1; break;	}
			}
			if(!iFlag)
			{
				VerIdx vi; vi.Ver = *r; vi.idx = index; vi.times = times;
				oddpoints->push_back(vi);
			}
		}
		index++;
	}
}


#define THOLD 5
void TerrainProcessingUtils::computePerimeterAll( osg::Geometry *g, std::vector<osg::ref_ptr<osg::Vec3Array>> &perimeterList )
{
	std::vector<Edge>   edges;
	osg::TriangleFunctor<MakeEdgeList> mel;
	mel.setEdgeVector(&edges);
	g->accept( mel );
	
	/* �ж�edges�������еıߣ��������߳��ֶ�Σ��Ͱ�����߱�עΪfalse����ʾ����������ڲ��ıߡ�2012-10-10 */
	for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		std::vector<Edge>::iterator q = p;
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
	std::vector<Edge> kept;
	for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
	{
		if(( p->keep ) && (p->A != p->B ))
		kept.push_back( *p );
	}


	/* ��һ�γ���Ϊ���ҵ��߻�����Ķ�ͷ�㣬�����������к��ʵĶ�ͷ�㡣*/
	std::vector<VerIdx> OddPoints;
	std::vector<Edge> OddPointPairs;

	/* ��һ��������㣬�ж�ÿ������ֵĴ�����������ֵĴ�����������˵���ж�ͷ�㡣2012-10-10 */
	if(1)
	{
		/* ��һ��Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ�������1�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		FindOutAllOddPoints(g, &kept, &OddPoints, 0);

		/* */	
		if(OddPoints.size() > 1)
		{
			/* �ڶ�����Ҫ�ҳ����ԣ�������������㣬����������ͬ�������ֵ�λ�ò�ͬ��*/
			for(unsigned int i=0; i<OddPoints.size() - 1; i++)
			{
				for(unsigned int j=i+1; j<OddPoints.size(); j++)
				{
					if(OddPoints[i].Ver != OddPoints[j].Ver)
					{
						Edge e(OddPoints[i].Ver, OddPoints[j].Ver); OddPointPairs.push_back(e);	
					}	
				}
			}

			/* ���������������Դ�ԭ�б߼��������ҵ���Ӧ�ıߣ����뵽kept�߼������档 */
			if(OddPointPairs.size() > 0)
			{
				for(std::vector<Edge>::iterator q = OddPointPairs.begin(); q != OddPointPairs.end(); q++)
				{
					/* �ҳ���һ���߻� */
					int iFlag = 0;
					for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
					{
						if(((p->A == q->A ) && (p->B == q->B )) || ((p->A == q->B ) && (p->B == q->A )))
						{
							kept.push_back( *q );	q->keep = true;
							break;
						}
					}
				}
			}

			/* �����������ͷ�㣬��ô�������Ļ���Ϣ�ߣ������ǽ�������*/
			/* ��һ��Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ�������1�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
			OddPoints.clear();	OddPointPairs.clear();
			FindOutAllOddPoints(g, &kept, &OddPoints, 0);
			if(OddPoints.size() == 2)
			{
				if(perimeterList.size() > 0)
				{
					for( std::vector<osg::ref_ptr<osg::Vec3Array>>::iterator q= perimeterList.begin(); q != perimeterList.end(); q++ )
					{
						/* ���ȣ�������߻����ҵ�һ������� */
						for( osg::Vec3Array::iterator pvx = q->get()->begin(); pvx != q->get()->end(); pvx++)
						{
							osg::Vec3 qSpace = *pvx;
							
							int iFlag = 0, curOPI = 0;
							for(unsigned int i=0; i<OddPoints.size(); i++)
							{
								if(qSpace == OddPoints[i].Ver)
								{
									iFlag = 1; curOPI = i;	break;
								}
							}
							
							/* ���߻����ҵ����������� */
							if(iFlag == 1)
							{
								/* ���ȱ�֤��������һ���߶Σ�����kept�����ڡ�*/
								int iFlag2 = 0;
								osg::Vec3 qNext = *(pvx + 1);
								Edge e(qSpace, qNext);
								for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
								{
									if(((p->A == e.A ) && (p->B == e.B )) || ((p->A == e.B ) && (p->B == e.A )))
									{
										iFlag2 = 1;		break;
									}
								}
								if(iFlag2 == 1) continue;
								
								/* ����߶Σ�ֱ����������һ������� */
								iFlag2 = 0;
								osg::Vec3Array::iterator ir = pvx;
								for(; ir != q->get()->end(); ir++)
								{
									if(*ir == OddPoints[1-curOPI].Ver)
									{
										iFlag2 = 1; break;
									}
									
									osg::Vec3 qN = *(ir + 1);
									Edge e(*ir, qN);
									kept.push_back(e);
								} 
								if(iFlag2) break;
							}
						}
					}
				}
			}
		}else if(OddPoints.size() == 1)
		{

			/* �� �������� ����ߣ���ԭ�� kept �����õ���Ȼ�����¼��㡣 */
			while(OddPoints.size() > 0)
			{
				std::vector<Edge> remain;
				for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
				{
					if (!((p->A == OddPoints[0].Ver) || (p->B == OddPoints[0].Ver)))	remain.push_back(*p);
				}
				kept.clear();	OddPoints.clear();
				for( std::vector<Edge>::iterator p = remain.begin(); p != remain.end(); p++ )
				{
					p->keep = true;
					kept.push_back( *p );
				}
				
				/* ������㣬����ҵ��ͼ���ȥ��� */
				FindOutAllOddPoints(g, &kept, &OddPoints, 0);
			}

		}
	}


	/* �ڶ�������㣬�Ҷ���ʹ����Ϊ����2�������εġ�*/
	OddPoints.clear();	OddPointPairs.clear();
	if(1)
	{
		/* ��һ��Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		FindOutAllOddPoints(g, &kept, &OddPoints, 1);

		/* ���裺�������������֮�䣬�й��������ߣ���ô�����й�����������ȥ���� */
		if(OddPoints.size() > 1)
		{
			/* �ڶ�����Ҫ�ҳ����ԣ�������������㣬����������ͬ�������ֵ�λ�ò�ͬ��*/
			for(unsigned int i=0; i<OddPoints.size() - 1; i++)
			{
				for(unsigned int j=i+1; j<OddPoints.size(); j++)
				{
					if(OddPoints[i].Ver != OddPoints[j].Ver)
					{
						Edge e(OddPoints[i].Ver, OddPoints[j].Ver); OddPointPairs.push_back(e);	
					}	
				}
			}

			/* ���������������Դ�ԭ�б߼��������ҵ���Ӧ�ıߣ����뵽kept�߼������档 */
			if(OddPointPairs.size() > 0)
			{
				for(std::vector<Edge>::iterator q = OddPointPairs.begin(); q != OddPointPairs.end(); q++)
				{
					/* �ҳ���һ���߻� */
					int iFlag = 0;
					for( std::vector<Edge>::iterator p = edges.begin(); p != edges.end(); p++ )
					{
						if(((p->A == q->A ) && (p->B == q->B )) || ((p->A == q->B ) && (p->B == q->A )))
						{
							kept.push_back( *q );	q->keep = true;
							break;
						}
					}
				}
			}
		}
	}


	/* ����������㣬�Ҷ���ʹ����Ϊ����2��ż���εġ�*/
	OddPoints.clear();	OddPointPairs.clear();
#if 0	/* Deleted by Galen 2013-03-15 */
	if(1)
#else
	if(0)
#endif
	{
		/* ��һ��Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		FindOutAllOddPoints(g, &kept, &OddPoints, 1);

		/* ���裺�������������֮�䣬û�й��������ߡ� */
		if(OddPoints.size() > 0)
		{
			for(unsigned int i=0; i<OddPoints.size(); i++)
			{
				/* ��һ�����ռ�����ͬ�����ıߡ�*/
				std::vector<Edge> withOdd;
				for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
				{
					if ( ((p->A == OddPoints[i].Ver) || (p->B == OddPoints[i].Ver)))	withOdd.push_back(*p);
				}
				
				/* �ڶ���������ÿһ�����ߣ����ռ�һ�����߻������ʮ�����ھͳ�Ϊ�߻���������������ߡ�*/
				for( std::vector<Edge>::iterator pp = withOdd.begin(); pp != withOdd.end(); pp++ )
				{
					/* �ҵ����ߵ���㡣*/
					VerIdx currOddPnt;
					for(unsigned int j=0; j<OddPoints.size(); j++)
					{
						if((OddPoints[j].Ver == pp->A)	|| (OddPoints[j].Ver == pp->B))
						{	currOddPnt = OddPoints[j];	break;	}
					}
					
					/* ��������߿�ʼ�����߻� */
					Edge start = *pp;
					Edge e = start;
					osg::Vec3 A = currOddPnt.Ver;
					osg::Vec3 B; if(e.A == A) B = e.B; else B = e.A;
					osg::Vec3 firstA = A;
					int iCounter = 0;
					bool bEnd = true;	//�����ǰ�����бߴղ���һ���߻���bEnd�ͻ����true�����������˳�do..whileѭ��������ܴճ�һ���߻���bEndʼ��Ϊfalse.

					/* ����kept����� start ����ߣ��� false��*/
					for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
					{
						if(*p == start)	{ p->keep = false; break;	}
					}

					/* �ٽ��ż����� */
					do 
					{
						bEnd = true;
						for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
						{
							// �ڲ���ͨ�û���ʱ����������������ԭ��δ֪���п���������Ĵ����󵥱߳�����Ϊ������һ����������������Ĺ��ߣ�Ҳ�����ɵ���, ������Щ�ж���Ϊ�˱�������
							if( p->has(B) && !(e == *p) && p->keep )  
							{
								A = B;
								if( B == p->A )
									B = p->B;
								else
									B = p->A;
			
								p->keep = false;
								e = *p;
								bEnd = false;

								/* ���������ж��˳����� */
								iCounter++;
								if((iCounter >= THOLD) || (p->has(firstA)))
								{
									bEnd = true;	break;
								}
								
								break;
							}
						}
						if (bEnd)
						{
							if((start.has(e.A)) || (start.has(e.B)))
							{	
								e = start;	bEnd = false;
							}
						}
					}while( !(e == start) && !bEnd );

					/* �������С�� THOLD �Σ�������߻�ȥ����������һ���ߵļ��� */
					if(iCounter < THOLD)
					{
						std::vector<Edge> remain;
						std::vector<Edge> remove;	
						for( std::vector<Edge>::iterator q = kept.begin(); q != kept.end(); q++ )
						{
							if (q->keep == true)	remain.push_back(*q); else remove.push_back(*q);
						}
						kept.clear();
						for( std::vector<Edge>::iterator q = remain.begin(); q != remain.end(); q++ )
						{
							kept.push_back( *q );
						}

						/* ������е�kept��û�ж�ͷ�㣬����ж�ͷ�㣬��ôҪ��ɾ���ı�ȫ�ָ�������*/
						std::vector<VerIdx> Singles;
						FindOutAllOddPoints(g, &kept, &Singles, 0);
						if(Singles.size() > 0)
						{
							for( std::vector<Edge>::iterator q = remove.begin(); q != remove.end(); q++ )
							{
								kept.push_back( *q );
							}
						}
						
					}else{
						for( std::vector<Edge>::iterator q = kept.begin(); q != kept.end(); q++ )
							q->keep = true;
					}
					
					unsigned int kk = 0;
				}// for withOdd
			}// for OddPoints 
		}//if oddPoints;
	}// if 1

	/* �ѵ�ǰ�߼��ϱ������� */
	perimeterList.clear();
	std::vector<Edge> store;
	for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )		store.push_back( *p );
	int needLoop;
	std::vector<osg::ref_ptr<osg::Vec3Array>> localList;

CPA_LOOP:
	/* �ָ��߼��� */
	kept.clear();	needLoop = 0;
	for( std::vector<Edge>::iterator p = store.begin(); p != store.end(); p++ )		kept.push_back( *p );

	/* */
	while ( kept.size() > 2 )
	{
		osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;

CPA_LOOP_2:
		Edge start = kept.front();

		/* ���ҵ�����㣬���������㣬��ô������㿪ʼ�����û�У��ʹӵ�һ���㿪ʼ */
		OddPoints.clear();
		FindOutAllOddPoints(g, &kept, &OddPoints, 1);
		if(OddPoints.size() > 0)
		{
			for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
			{
				if (p->A == OddPoints[0].Ver)
				{
					start = *p;	break;
				}
			}
		}
		
		/* ���� */
		Edge e = start;
		osg::Vec3 A = e.A;
		osg::Vec3 B = e.B;
		osg::Vec3 firstA = A;
		bool bEnd = true;	//�����ǰ�����бߴղ���һ���߻���bEnd�ͻ����true�����������˳�do..whileѭ��������ܴճ�һ���߻���bEndʼ��Ϊfalse.
		do {
			bEnd = true;
			for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
			{
				// �ڲ���ͨ�û���ʱ����������������ԭ��δ֪���п���������Ĵ����󵥱߳�����Ϊ������һ����������������Ĺ��ߣ�Ҳ�����ɵ���, ������Щ�ж���Ϊ�˱�������
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
		if(e == start)
		{
			coords->push_back(firstA);
		}
		else
		{
			/* �� e ����ߣ���ԭ�� kept �����õ���Ȼ�����¼��㡣 */
			std::vector<Edge> remain;
			for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
			{
				if (!(*p == e))	remain.push_back(*p);
			}
			kept.clear();
			for( std::vector<Edge>::iterator p = remain.begin(); p != remain.end(); p++ )
			{
				p->keep = true;
				kept.push_back( *p );
			}

			/* ������LINELOOP */
			coords->clear();
			goto CPA_LOOP_2;
		}
		
		/* ѹ���߻����� */
		if(coords->size() >= 3)	perimeterList.push_back(coords);

		/* ���������Ϊ1��˵��ֻѹ����firstA�㣬Ҳ���� start�� ��һ����*/
		if(coords->size() == 2)
		{
			kept.front().keep = false;
		}

		/* ����Ƿ���ʣ��ı�û������µ��߻���*/
		std::vector<Edge> remain;
		for( std::vector<Edge>::iterator p = kept.begin(); p != kept.end(); p++ )
		{
			if (p->keep)
				remain.push_back(*p);
		}

		kept.clear();
		for( std::vector<Edge>::iterator p = remain.begin(); p != remain.end(); p++ )
		{
			kept.push_back( *p );
		}
	}// while(kept)


	/* Add by Galen 2013-02-22 �Ȱ�ֻ��һ���������ε��߻����С�*/
	if((perimeterList.size() == 1) && (perimeterList[0]->size() == 4))
		return;
	/* Add by Galen 2013-02-22 End */

	/* �ж������߻����Ƿ��зǷ��߻� 2012-10-12 */
	/* ����Ƿ��߻�����������"��Ȼ��һ���߻�����ʵ����ֻ��һ����"���߻���������ߵĳ��ȵ��ڻ�ӽ�����������߳���֮�͡�*/
	for(unsigned int i=0; i<perimeterList.size(); i++)
	{
		/* ���ж�����߻��Ƿ�Ϸ� */
		Edge* curr = CheckLINELOOPValibility(perimeterList[i]);
		if(curr != NULL)
		{
			/* �� kept ����ɾ�����"��Ȼ��һ���߻�����ʵ����ֻ��һ����"�������߻� */
			std::vector<Edge> remain;

			/* ��һ�����ҳ�����߻������б� */
			std::vector<Edge>   peri_edges;
			for( osg::Vec3Array::iterator p = perimeterList[i]->begin(); p != perimeterList[i]->end(); p++ )
			{
				if( p == perimeterList[i]->end() - 1 ) break;
				osg::Vec3Array::iterator q = p + 1;
				Edge e(*p, *q);	peri_edges.push_back(e);
			}

			/* �ڶ���������kept��������ȥ�����������ҵ��ı� */
			for( std::vector<Edge>::iterator p = store.begin(); p != store.end(); p++ )
			{
				int iFlag = 0;
				for( std::vector<Edge>::iterator q = peri_edges.begin(); q != peri_edges.end(); q++ )
				{
					if( *p == *q )	{ iFlag = 1;	break;	}
				}
				if (!iFlag)	remain.push_back(*p);
			}
			store.clear();
			for( std::vector<Edge>::iterator p = remain.begin(); p != remain.end(); p++ )
			{
				p->keep = true;
				store.push_back( *p );
			}

			needLoop = 1;
		}else{
			localList.push_back(perimeterList[i]);
		}
	}
	if(needLoop)
	{
		perimeterList.clear();
		goto CPA_LOOP;
	}
	
	/* �˳�����,���ʲôҲûȡ������ʹ��ԭ��ǰ���ֻ�õ�ֵ */
	if(perimeterList.size() == 0)
	{
		for(unsigned int i=0; i<localList.size(); i++)
			perimeterList.push_back(localList[i]);
	}

}

TerrainProcessingUtils::Edge::Edge(const Edge  &e ): 
    A(e.A), B(e.B), keep(true) 
{}

TerrainProcessingUtils::Edge::Edge(const osg::Vec3 &a, const osg::Vec3 &b ): 
    A(a), B(b) 
{}
        
bool TerrainProcessingUtils::Edge::operator == (Edge &rhs )
{
    if( (A == rhs.A && B == rhs.B) || (A == rhs.B && B == rhs.A))
        return true;
    else
        return false;
}
        
bool TerrainProcessingUtils::Edge::has( osg::Vec3 v )
{
    if( v == A || v == B ) return true;
    return false;
}

double TerrainProcessingUtils::Edge::length()
{
    return sqrt(((A.x()-B.x()) * (A.x()-B.x())) + ((A.y()-B.y()) * (A.y()-B.y())) + ((A.z()-B.z()) * (A.z()-B.z())));
}

TerrainProcessingUtils::MakeEdgeList::MakeEdgeList():
    edges(0L) {}
    
void TerrainProcessingUtils::MakeEdgeList::setEdgeVector( std::vector<Edge> *e )
{
    edges = e;
}
        
void TerrainProcessingUtils::MakeEdgeList::operator() (const osg::Vec3& v1,const osg::Vec3& v2,const osg::Vec3& v3, bool flag) const
{
    if( edges != 0L )
    {
        edges->push_back(Edge( v1, v2 ));
        edges->push_back(Edge( v2, v3 ));
        edges->push_back(Edge( v3, v1 ));
    }
}

