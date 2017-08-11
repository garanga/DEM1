/*************************************************************/
/*                                                           */
/* This program code belongs to the book                     */
/* "Computational Granular Dynamics" by                      */
/* Thorsten Poeschel and Thomas Schwager.                    */
/* ISBN: 3-540-21485-2                                       */
/* http://www.springeronline.com/3-540-21485-2               */
/*                                                           */
/* If you use the code for scientific purposes, please,      */
/* cite it as                                                */
/* "Computational Granular Dynamics"                         */
/* by T. Poeschel and T. Schwager.                           */
/* Springer (Berlin, New York, 2005). 322pp.                 */
/* For commercial use, please contact the authors.           */
/*                                                           */
/* Copyright (C) 2003,                                       */
/* Thorsten Poeschel (thorsten.poeschel@charite.de)          */
/* Thomas Schwager (thomas.schwager@charite.de)              */
/* This program code is distributed without any warranty.    */
/*                                                           */
/*************************************************************/

/*************************************************************/
/*                                                           */
/* This program code belongs to the book                     */
/* "Computational Granular Dynamics" by                      */
/* Thorsten Poeschel and Thomas Schwager.                    */
/* ISBN: 3-540-21485-2                                       */
/* http://www.springeronline.com/3-540-21485-2               */
/*                                                           */
/* If you use the code for scientific purposes, please,      */
/* cite it as                                                */
/* "Computational Granular Dynamics"                         */
/* by T. Poeschel and T. Schwager.                           */
/* Springer (Berlin, New York, 2005). 322pp.                 */
/* For commercial use, please contact the authors.           */
/*                                                           */
/* Copyright (C) 2003,                                       */
/* Thorsten Poeschel (thorsten.poeschel@charite.de)          */
/* Thomas Schwager (thomas.schwager@charite.de)              */
/* This program code is distributed without any warranty.    */
/*                                                           */
/*************************************************************/

#ifndef _Sphere_h
#define _Sphere_h

#include <iostream>
#include "Vector.h"
using namespace std;

inline double
normalize(double dx, double L)
{
    while(dx<-L/2) dx+=L;
    while(dx>=L/2) dx-=L;
    return dx;
}

class Sphere
{
    friend istream&
    operator >> (istream&, Sphere&);

    friend ostream&
    operator << (ostream&, const Sphere&);

    friend double
    Distance(const Sphere& p1, const Sphere& p2, double lx, double ly)
    {
        double dx = normalize(p1.rtd0_.x()-p2.rtd0_.x(),lx);
        double dy = normalize(p1.rtd0_.y()-p2.rtd0_.y(),ly);
        return sqrt(dx*dx+dy*dy);
    }

    friend void
    force(Sphere&, Sphere&, double, double);

public:

    Sphere():
        rtd0_(null), rtd1_(null), rtd2_(null), rtd3_(null), rtd4_(null)
    {

    }

    Vector&
    pos()
    {
        return rtd0_;
    }

    Vector
    pos() const
    {
        return rtd0_;
    }

    double&
    x()
    {
        return rtd0_.x();
    }

    double
    x() const
    {
        return rtd0_.x();
    }

    double&
    y()
    {
        return rtd0_.y();
    }

    double
    y() const
    {
        return rtd0_.y();
    }

    double&
    phi()
    {
        return rtd0_.phi();
    }

    double
    phi() const
    {
        return rtd0_.phi();
    }

    double&
    vx()
    {
        return rtd1_.x();
    }

    double
    vx() const
    {
        return rtd1_.x();
    }

    double&
    vy()
    {
        return rtd1_.y();
    }

    double
    vy() const
    {
        return rtd1_.y();
    }

    double&
    omega()
    {
        return rtd1_.phi();
    }

    double
    omega() const
    {
        return rtd1_.phi();
    }

    const Vector&
    velocity() const
    {
        return rtd1_;
    }

    double&
    r()
    {
        return r_;
    }

    double
    r() const
    {
        return r_;
    }

    double
    m() const
    {
        return m_;
    }

    int
    ptype() const
    {
        return ptype_;
    }

    void
    predict(double);

    void
    add_force(const Vector& f)
    {
        force_+=f;
    }

    void
    correct(double);

    double
    kinetic_energy() const;

    void
    set_force_to_zero()
    {
        force_ = null;
    }

    void
    boundary_conditions(int, double, double);

    void
    periodic_bc(double, double, double, double);

private:
    double m_, r_, J_;
    int ptype_;
    double Y_, mu_, A_, gamma_;
    Vector rtd0_, rtd1_, rtd2_, rtd3_, rtd4_;
    Vector force_;
};

#endif


