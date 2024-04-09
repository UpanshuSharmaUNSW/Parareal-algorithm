#include"particle.hpp"

particle::particle(const input & I)
{
    // // ** Double well: initial configuration is bottom of left well with some initial momentum
    // q = -1.;
    // p = 0.1;

    // ** Single well: minimum of well
    q = 1.;
    p = 0.;

    // // ** Lennard-Jones: minimum of LJ well is at x=1 (LJ(1)=-1)
    // q = 1.05;
    // p = 0.;
}

particle & particle::operator= (const particle & a)
// affectation surchargée
{
    if (this!=&a)
    {
        //	cerr<<"faire qqch dans ="<<endl;
        q = a.q;
        p = a.p;
    }
    else
    {
        //	cerr<<"ne rien faire dans ="<<endl;
    }
    return *this;
}
