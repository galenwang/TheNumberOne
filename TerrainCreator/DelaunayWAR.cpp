//#define DEBUG 
#ifdef DEBUG
#include <osg/Geode>
#include <osg/Geometry>
#include <osgDB/WriteFile>
#endif

#include <osg/Notify>
#include "TerrainProcessingUtils.h"
#include "DelaunayWAR.h"

#include <osg/Vec2>

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

DelaunayWAR::DelaunayWAR( osg::Vec3Array *coords) 
{
    _coords = coords;
}

bool DelaunayWAR::triangulate()
{
    bool flag = true;
    const unsigned int MAX_TRIES = 3;
    unsigned int tries = 0;
    do {
        _dt = new osgUtil::DelaunayTriangulator(_coords.get());
        _dt->triangulate();
        flag = _apply();
        tries++;
        if( tries > MAX_TRIES )
        {
            osg::notify(osg::WARN) << "Terrain Deformation Software:  DelaunayWAR::triangulate() - MAX_TRIES reached" << std::endl;
            flag = false;
        }
    } while( flag );
    return flag;
}

const osg::DrawElementsUInt *DelaunayWAR::getTriangles() const
{
    if( !_dt.valid() )
        return 0L;
    return _dt->getTriangles();
}

osg::DrawElementsUInt *DelaunayWAR::getTriangles()
{
    if( !_dt.valid() )
        return 0L;
    return _dt->getTriangles();
}

//bool DelaunayWAR::_apply( osg::Vec3Array &coords, osgUtil::DelaunayTriangulator &dt )
bool DelaunayWAR::_apply()
{
    // For every pair of points on the periphery, find at least one triangle that has
    // both points.  Record any pair of points that are don't have at least one 
    // triangle with both points
    //osg::ref_ptr<osg::Vec3Array> new_coords = _periphery( *(_coords.get()) );

    osg::ref_ptr<osg::Vec3Array> new_coords = TerrainProcessingUtils::computePeriphery( *(_coords.get()) );

#ifdef DEBUG
    { // DeBUG - view periphery

        osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
        colors->push_back( osg::Vec4( 0.0, 0.0, 1.0, 1.0 ));

        osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;
        geometry->setVertexArray( new_coords.get() );
        geometry->setColorBinding( osg::Geometry::BIND_OVERALL );
        geometry->setColorArray( colors.get() );
        geometry->addPrimitiveSet( new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, new_coords->size()));

        osg::ref_ptr<osg::Geode> geode = new osg::Geode;
        geode->addDrawable( geometry.get() );
		osgDB::writeNodeFile( *geode.get(), "E:\\delWarPeriph.osg" );
    }
#endif

    std::vector<PairPts> nonpairs;
    unsigned int n = new_coords->size();
    for( unsigned int i = 0; i < n; i++ ) 
    {
        osg::Vec3 A = (*new_coords.get())[i];
        osg::Vec3 B = (*new_coords.get())[(i+1)%n];

        if( !_pairInTriangles( A, B, *(_coords.get()), _dt->getTriangles() ) )
            nonpairs.push_back( PairPts( A, B ) );
    }

#ifdef DEBUG
    printf("DEBUG: sizeof nonpairs: %d\n", nonpairs.size() );
    for( std::vector<PairPts>::iterator p = nonpairs.begin();
            p != nonpairs.end(); p++ )
    {
        printf( "         %f %f %f  -  %f %f %f\n", p->A[0], p->A[1], p->A[2], p->B[0], p->B[1], p->B[2] );
    }
#endif

    if( nonpairs.size() > 0 )
    {
        for( std::vector<PairPts>::iterator ip = nonpairs.begin(); 
             ip != nonpairs.end(); 
             ip++ )
        {
            osg::Vec3 C = ((*ip).A + (*ip).B) * 0.5;
#ifdef DEBUG
            printf("DEBUG: adding a vector (%f %f %f)\n", C[0], C[1], C[1] );
            printf("DEBUG: Before adding sizeof coords is %d\n", _coords->size() );
#endif
            _coords->push_back( C );
#ifdef DEBUG
            printf("DEBUG: After adding sizeof coords is %d\n", _coords->size() );
#endif
        }
        return true;
    }
    else
        return false;

#if 0
    // 2. For each point in the pairs: find all triangles they belong to
    //    and list all points in each triangle which is not "this" point.
    //    This should result in two lists.  Compare both lists for common points
    osg::ref_ptr<osg::Vec3Array> addTriangles = new osg::Vec3Array;
    for( std::vector<PairPts>::iterator ip = nonpairs.begin(); 
             ip != nonpairs.end(); 
             ip++ )
    {
        osg::Vec3 A = ip->A;
        osg::Vec3 B = ip->B;
        std::set<osg::Vec3> listA = 
            _vec3ListOfTrianglesPtBelongsTo(A, coords, dt.getTriangles() );
        std::set<osg::Vec3> listB = 
            _vec3ListOfTrianglesPtBelongsTo(B, coords, dt.getTriangles() );

        bool found = false;
        osg::Vec3 C;
        for( std::set<osg::Vec3>::iterator pa = listA.begin(); pa != listA.end(); pa++ )
        {
            for( std::set<osg::Vec3>::iterator pb = listB.begin(); pb != listB.end(); pb++ )
            {
                if( (*pa) == (*pb) )
                {
                    C = *pa;
                    found = true;
                    break;
                }
            }
            if( found )
                break;
        }

        if( found )
        {
            addTriangles->push_back( A ); 
            addTriangles->push_back( B ); 
            addTriangles->push_back( C ); 
        }
    }

    for( osg::Vec3Array::iterator p = addTriangles->begin();
            p != addTriangles->end();
            p++  )
    {
        unsigned int index = coords.size();
        coords.push_back( *p );
        dt.getTriangles()->push_back( index );
    }
#endif
}

osg::ref_ptr<osg::Vec3Array> DelaunayWAR::_periphery( const osg::Vec3Array &v )
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
    osg::Vec2 C( (minX + maxX)/2, (minY + maxY)/2 );

    //  Find the point furthest (F) from the center
    osg::Vec3 F;
    double dist = 0.0;
    for( p = v.begin(); p != v.end(); p++ )
    {
        osg::Vec3 v3 = *p;
        osg::Vec2 v2(v3[0], v3[1]);
        //double l = (*p - C).length();
        double l = (v2 - C).length();
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
    osg::Vec3 L = osg::Vec3(C[0],C[1], 0);
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
            if( bb  <= smallestAngle )
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


bool DelaunayWAR::_pairInTriangles( 
        const osg::Vec3 &A, const osg::Vec3 &B,  
        const osg::Vec3Array &coords, 
        osg::DrawElementsUInt *t )
{
    osg::DrawElementsUInt::iterator p;
    bool found = false;
    for( p = t->begin() ; p != t->end(); p+=3 )
    {
        osg::Vec3 a = coords[*(p+0)];
        osg::Vec3 b = coords[*(p+1)];
        osg::Vec3 c = coords[*(p+2)];

        if( (A == a || A == b || A == c ) &&
            (B == a || B == b || B == c ) )
        {
            found = true;
            break;
        }
    }
    return found;
}

#if 0
std::set<osg::Vec3> DelaunayWAR::_vec3ListOfTrianglesPtBelongsTo( osg::Vec3 P, 
        const osg::Vec3Array &coords, 
        osg::DrawElementsUInt *t )
{
    std::set<osg::Vec3> list;
    for( osg::DrawElementsUInt::iterator p = t->begin() ; p != t->end(); p+=3 )
    {
        osg::Vec3 a = coords[*(p+0)];
        osg::Vec3 b = coords[*(p+1)];
        osg::Vec3 c = coords[*(p+2)];

        if( P == a )
        {
            list.insert( b );
            list.insert( c );
        }
        else if( P == b )
        {
            list.insert( a );
            list.insert( c );
        }
        else if( P == c )
        {
            list.insert( a );
            list.insert( b );
        }
    }

    return list;
}
#endif



