#include "fractals.hpp"

void init_fractal(const size_t& x, const size_t& y, const size_t& width, const size_t& height, const size_t& square_scale, complex<long double>& z, complex<long double>& c, [[maybe_unused]] vector<complex<long double>>& history, [[maybe_unused]] vector<complex<long double>>& extra_params, const fractalsettings_t& fset){
    switch(fset.fractal_type){
        default:
        case ftype::mandelbrot:
        case ftype::tippets:
        case ftype::burnship:
        case ftype::mandelbar:
        case ftype::magnet1:
        case ftype::magnet2:
        case ftype::cactus0:
        case ftype::unpolsquare:
        case ftype::moth:
        case ftype::cave:
        case ftype::wankel:
        case ftype::seaangel:
        case ftype::smith:

        case ftype::spade:
            if(fset.julia_mode){
                z.real(real_from_posX(x, width, square_scale, fset));
                z.imag(imag_from_posY(y, height, square_scale, fset));
                c.real(fset.julia_re);
                c.imag(fset.julia_im);
            } else {
                z.real(0);
                z.imag(0);
                c.real(real_from_posX(x, width, square_scale, fset));
                c.imag(imag_from_posY(y, height, square_scale, fset));
            }
            break;
            
        case ftype::cactus:
        case ftype::zubieta:
        case ftype::zubitheta:
            if(fset.julia_mode){
                z.real(real_from_posX(x, width, square_scale, fset));
                z.imag(imag_from_posY(y, height, square_scale, fset));
                c.real(fset.julia_re);
                c.imag(fset.julia_im);
            } else {
                z.real(real_from_posX(x, width, square_scale, fset));
                z.imag(imag_from_posY(y, height, square_scale, fset));
                c.real(real_from_posX(x, width, square_scale, fset));
                c.imag(imag_from_posY(y, height, square_scale, fset));
            }
            break;

        case ftype::logmap:
            if(fset.julia_mode){
                z.real(real_from_posX(x, width, square_scale, fset));
                z.imag(imag_from_posY(y, height, square_scale, fset));
                c.real(fset.julia_re);
                c.imag(fset.julia_im);
            } else {
                z.real(0.5l);
                z.imag(0);
                c.real(real_from_posX(x, width, square_scale, fset));
                c.imag(imag_from_posY(y, height, square_scale, fset));
            }
            break;
    }

    return;
}