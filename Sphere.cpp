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

#include "Sphere.h"
#include <assert.h>

extern Vector G;

void
force(Sphere& p1, Sphere& p2, double lx, double ly)
{
    double dx=normalize(p1.x()-p2.x(),lx);
    double dy=normalize(p1.y()-p2.y(),ly);
    double rr=sqrt(dx*dx+dy*dy);
    double r1=p1.r();
    double r2=p2.r();
    double xi=r1+r2-rr;

    if(xi>0) // The particles overlap
    {
        double Y=p1.Y_*p2.Y_/(p1.Y_+p2.Y_); // nu is zero
        double A=0.5*(p1.A_+p2.A_);
        double mu = (p1.mu_<p2.mu_ ? p1.mu_ : p2.mu_);
        double gamma = (p1.gamma_<p2.gamma_ ? p1.gamma_ : p2.gamma_);
        double reff = (r1*r2)/(r1+r2);
        double dvx=p1.vx()-p2.vx();
        double dvy=p1.vy()-p2.vy();
        double rr_rez=1/rr;
        double ex=dx*rr_rez;
        double ey=dy*rr_rez;
        double xidot=-(ex*dvx+ey*dvy);
        double vtrel=-dvx*ey + dvy*ex + p1.omega()*p1.r()-p2.omega()*p2.r();
        double fn=sqrt(xi)*Y*sqrt(reff)*(xi+A*xidot);
        double ft=-gamma*vtrel;

        if(fn<0)      fn=0;
        if(ft<-mu*fn) ft=-mu*fn;
        if(ft>mu*fn) ft=mu*fn;
        if(p1.ptype()==0)
        {
            p1.add_force(Vector(fn*ex-ft*ey, fn*ey+ft*ex, r1*ft));
        }
        if(p2.ptype()==0)
        {
            p2.add_force(Vector(-fn*ex+ft*ey, -fn*ey-ft*ex, -r2*ft));
        }
    }
}

void
Sphere::predict(double dt)
{
    double a1=dt;
    double a2=a1*dt/2;
    double a3=a2*dt/3;
    double a4=a3*dt/4;

    rtd0_ += a1*rtd1_ + a2*rtd2_ + a3*rtd3_ + a4*rtd4_;
    rtd1_ += a1*rtd2_ + a2*rtd3_ + a3*rtd4_;
    rtd2_ += a1*rtd3_ + a2*rtd4_;
    rtd3_ += a1*rtd4_;
}

void
Sphere::correct(double dt)
{
    static Vector accel,corr;
    double dtrez = 1.0/dt;
    const double coeff0 = 19.0/90.0*( dt*dt/2.0);
    const double coeff1 =  3.0/ 4.0*(    dt/2.0);
    const double coeff3 =  1.0/ 2.0*( 3.0*dtrez);
    const double coeff4 =  1.0/12.0*(12.0*(dtrez*dtrez));

    accel = Vector(
            (1.0/m_)*force_.x()   + G.x(),
            (1.0/m_)*force_.y()   + G.y(),
            (1.0/J_)*force_.phi() + G.phi());
    corr=accel-rtd2_;
    rtd0_ += coeff0*corr;
    rtd1_ += coeff1*corr;
    rtd2_ = accel;
    rtd3_ += coeff3*corr;
    rtd4_ += coeff4*corr;
}

void
Sphere::boundary_conditions(int n, double timestep, double Time)
{
    switch(ptype())
    {
        case(0): break;     // normal granular particle
        case(1): break;     // static wall
        case(2):            // oscillating wall
        {
            x() = 0.5-0.4*cos(10.0*Time);
            y() = 0.1;
            vx() = 10.0*0.4*sin(Time);
            vy() = 0.0;
        } break;
        case(3):            // rotating wall
        {
            double xx = x()-0.5;
            double yy = y()-0.5;
            double xp = xx*cos(timestep) - yy*sin(timestep);
            double yp = xx*sin(timestep) + yy*cos(timestep);

            x() = 0.5 + xp;
            y() = 0.5 + yp;
            vx() = -yp;
            vy() =  xp;
            omega()=1;
        } break;

        case(4):
        {
            x() = 0.5+0.1*cos(Time) + 0.4*cos(Time+2*n*M_PI/128);
            y() = 0.5+0.1*sin(Time) + 0.4*sin(Time+2*n*M_PI/128);
            vx() =-0.1*sin(Time) - 0.4*sin(Time+2*n*M_PI/128);
            vy() = 0.1*cos(Time) - 0.4*cos(Time+2*n*M_PI/128);
            omega()=1;
        } break;

        case(5):
        {
            y() = 0.1+0.02*sin(30*Time);
            vx() = 0;
            vy() = 0.02*30*cos(30*Time);
        } break;

        case(6):
        {
            int i = n/2;
            y() = i*0.02+0.1+0.02*sin(30*Time);
            vx() = 0;
            vy() = 0.02*30*cos(30*Time);
        } break;

        default:
        {
            cerr << "ptype: " << ptype() << " not implemented\n";
            abort();
        }
    }
}

void
Sphere::periodic_bc(double x_0, double y_0, double lx, double ly)
{
    while(rtd0_.x()<x_0   ) rtd0_.x()+=lx;
    while(rtd0_.x()>x_0+lx) rtd0_.x()-=lx;
    while(rtd0_.y()<y_0   ) rtd0_.y()+=ly;
    while(rtd0_.y()>y_0+ly) rtd0_.y()-=ly;
}

double
Sphere::kinetic_energy() const
{
    return m_*(rtd1_.x()*rtd1_.x()/2.0 + rtd1_.y()*rtd1_.y()/2.0) +
           J_*rtd1_.phi()*rtd1_.phi()/2.0;
}

istream&
operator >> (istream& is, Sphere& p)
{
  is >> p.rtd0_ >> p.rtd1_
     >> p.r_ >> p.m_ >> p.ptype_
     >> p.Y_ >> p.A_ >> p.mu_ >> p.gamma_
     >> p.force_
     >> p.rtd2_ >> p.rtd3_ >> p.rtd4_;
  p.J_= p.m_*p.r_*p.r_/2.0;
  return is;
}

ostream&
operator << (ostream& os, const Sphere& p)
{
  os << p.rtd0_ << " " << p.rtd1_ << " ";
  os << p.r_ << " " << p.m_ << " " << p.ptype_ << " ";
  os << p.Y_ << " " << p.A_ << " " << p.mu_ << " " << p.gamma_ << " ";
  os << p.force_ << " ";
  os << p.rtd2_ << " " << p.rtd3_ << " " << p.rtd4_ << "\n" << flush;
  return os;
}
