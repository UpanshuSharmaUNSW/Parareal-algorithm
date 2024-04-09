#include"particle.hpp"

particle::particle()
{
    q_x = 1.0;
    q_y = 0.0;
    p_x = 0.0;
    p_y = 0.0;
}

particle::particle(double my_q_x, double my_q_y, double my_p_x, double my_p_y)
{
    q_x = my_q_x;
    q_y = my_q_y;
    p_x = my_p_x;
    p_y = my_p_y;
}

particle & particle::operator= (const particle & a)
// overload assignment
{
    if (this!=&a)
    {
        q_x = a.q_x;
        q_y = a.q_y;
        p_x = a.p_x;
        p_y = a.p_y;
    }
    else
    {
        //    cerr<<"do nothing ="<<endl;
    }
    return *this;
}
