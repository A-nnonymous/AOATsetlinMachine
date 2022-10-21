#include "Habitat.h"
#include <climits>

struct utility
{
    double value;
    utility()
    {
        value = -(__DBL_MAX__);
    }
    utility(double x)
    {
        value = x;
    }
};
struct point
{
    double x;
    double y;
    point()
    {
        x=0;y=0;
    }
    point(double xin, double yin)
    {
        x = xin;
        y = yin;
    }
};

utility concave(point p)
{
    double num = p.x*p.x + p.y*p.y;
    auto result = utility(-num);
    return result;
}


int main(int argc, char const *argv[])
{

    std::vector<double> min{-10,-10};
    std::vector<double> max{10,10};
    auto limits = Predator<utility,point,double>::rangeLimits(min,max);
    auto search = Predator<utility,point,double>::searchArgs(2,100,0.1,0.005,limits);
    auto test = Habitat<utility,point,double>(3,concave,point(), search);
    test.summonAll();
    return 0;
}
