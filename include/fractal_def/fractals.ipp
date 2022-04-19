#include "fractals.hpp"

inline long double real_from_posX(const size_t posX, const size_t width, const size_t squareScale, const fractalsettings_t fset){
    return 4.0l * ((long double)posX - (long double)width / 2.0l)  / ((long double)squareScale * fset.scaling_factor) + fset.offset_re;
}

inline long double imag_from_posY(const size_t posY, const size_t height, const size_t squareScale, const fractalsettings_t fset){
    return 4.0l * ((long double)posY - (long double)height / 2.0l) / ((long double)squareScale * fset.scaling_factor) + fset.offset_im;
}

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

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//MANDELBROT SET
//z(n+1) = z(n)^2 +c
template<typename T>
complex<T> mandelbrot(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    //z = re + i*im
    //c = cre + i*cim
    //z^2 + c = re^2 - im^2 + cre + i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = 2*re*im + cim

    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TIPPETS MANDELBROT SET
//z(n+1) = z(n)^2 +c but wrongly implemented
template<typename T>
complex<T> tippets_mandelbrot(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    //Normal:
    //im = 2 * re * im + cim;
    //re = re_2 - im_2 + cre;
    //
    //Tippets:
    //re = re_2 - im_2 + cre;
    //im = 2 * re * im + cim;

    const long double new_re = last_z.real() * last_z.real() - last_z.imag() * last_z.imag() + c.real();
    const long double new_im = 2*new_re*last_z.imag() + c.imag();

    const complex<T> new_z{new_re, new_im};
    /*
    complex<T> new_z = last_z * last_z + c;

    new_z.imag(new_z.imag() - 2 * (last_z.real() - last_z.real() * last_z.real() + last_z.imag() * last_z.imag() - c.real()) * last_z.imag());
    */

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//BURNING SHIP
//z(n+1) = (abs(re(z(n)) + i*abs(im(z(n)))) ^ 2 + c
template<typename T>
complex<T> burning_ship(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    //z = re + i*im
    //c = cre + i*cim
    //(abs(re) + i*abs(im))^2 + c = re^2 - im^2 + cre + abs(i*2*re*im) + i*cim;
    //re = re^2 - im^2 + cre
    //im = abs(2*re*im) + cim

    const complex<T> new_z = last_z * last_z + complex<T>{c.real(), (last_z.real() * last_z.imag() < 0 ? (-4 * last_z.real() * last_z.imag()) : (0)) + c.imag()};

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//MANDELBAR SET
//z(n+1) = conj(z(n))^2 +c
template<typename T>
complex<T> mandelbar(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    //z = re + i*im
    //c = cre + i*cim
    //conj(z)^2 + c = re^2 - im^2 + cre - i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = -2*re*im + cim

    const complex<T> new_z = conj(last_z) * conj(last_z) + c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//MAGNET FRACTAL TYPE I
//z(n+1) = ((z(n)^2 + c - 1) / (2*z(n) + c - 2))^2
template<typename T>
complex<T> magnet_type1(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> numerator = last_z * last_z + c - complex<T>(1, 0);
    const complex<T> denominator = complex<T>(2, 0) * last_z + c - complex<T>(2, 0);
    
    const complex<T> new_z = (numerator * numerator) / (denominator * denominator);

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//MAGNET FRACTAL TYPE II
//z(n+1) = ((z(n)^3 + 3*(c - 1) * z(n) + (c - 1)*(c - 2)) / (3z(n)^2 + 3*(c - 2)*z(n) + (c - 1)*(c - 2) + 1))^2
template<typename T>
complex<T> magnet_type2(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> numerator = last_z*last_z*last_z + complex<T>(3, 0) * (c - complex<T>(1, 0)) * last_z + (c - complex<T>(1, 0)) * (c - complex<T>(2, 0));
    const complex<T> denominator = complex<T>(3, 0) * last_z*last_z + complex<T>(3, 0) * (c - complex<T>(2, 0)) * last_z + (c - complex<T>(1, 0)) * (c - complex<T>(2, 0)) + complex<T>(1, 0);
    
    const complex<T> new_z = (numerator * numerator) / (denominator * denominator);

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//CACTUS FRACTAL
//CACTUS0 FRACTAL
//z(n+1) = z(n)^3 + (c - 1) * z(n) - c
template<typename T>
complex<T> cactus(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z*last_z*last_z + (c - complex<T>(1, 0)) * last_z - c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//ZUBIETA FRACTAL
//I don't know the real name of this fractal, I stole it from here http://paulbourke.net/fractals/Zubieta/
//z(n+1) = z(n)^2 + c / z(n)
template<typename T>
complex<T> zubieta(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c / last_z;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//ZUBITHETA FRACTAL
//Because after eta comes theta. This is very similar to the Zubieta fractal
//z(n+1) = z(n)^2 + z(n) / c
template<typename T>
complex<T> zubitheta(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + last_z / c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//LOGISTIC MAP FRACTAL
//z(n+1) = c * z(n) * (1-z(n))
template<typename T>
complex<T> logistic_map(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = c * last_z * (complex<T>(1, 0) - last_z);

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//UNPOLISHED SQUARE FRACTAL
//Re(n) = (Re(n-1) - Im(n-1) )*|Im(n-1)| + Re_c
//Im(n) = (Re(n-1) + Im(n-1) )*|Re(n-1)| + Im_c
template<typename T>
complex<T> unpol_square(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z{   (last_z.real() - last_z.imag()) * abs(last_z.imag()) + c.real(),
                                        (last_z.real() + last_z.imag()) * abs(last_z.real()) + c.imag()};

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//MOTH FRACTAL
//z(n+1) = (3*z^3 - 2*z - 1)/(z * c + 1)
template<typename T>
complex<T> moth(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z =
        (complex<T>(3,0) * last_z*last_z*last_z - complex<T>(2, 0) * last_z - complex<T>(1, 0)) / 
        (complex<T>(2, 0) * last_z * c + complex<T>(1, 0));

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//CAVE FRACTAL
//z(n+1) = (z(n)^3 + c)/(-2*z(n) + 1)
template<typename T>
complex<T> cave(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = (last_z*last_z*last_z + c) / (complex<T>(-2, 0)*last_z + complex<T>(1, 0));

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//WANKEL FRACTAL
//z(n+1) = (z(n)^3 + 1)/c
template<typename T>
complex<T> wankel(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = (last_z * last_z * last_z + complex<T>(1, 0)) / c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//SEA ANGEL FRACTAL
//z(n+1) = (z(n)^4 + 3*z(n)^2 + c) / (5*z(n)^2 - 3*z(n) + 2)
template<typename T>
complex<T> sea_angel(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z =
        (last_z*last_z*last_z*last_z + complex<T>(3,0)*last_z*last_z + c) /
        (complex<T>(5,0)*last_z*last_z - complex<T>(3,0)*last_z + complex<T>(2,0));

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//SMITH CHART FRACTAL
//z(n+1) = (1+z(n)) / (1-z(n)) + c
template<typename T>
complex<T> smith(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = (complex<T>(1, 0) + last_z) / (complex<T>(1, 0) - last_z) + c;

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//SPADE FRACTAL
//z(n+1) = z(n)^z(n) + z(n)/c
template<typename T>
complex<T> spadefract(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = pow(last_z, last_z) + last_z / c;    

    return new_z;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TEST FRACTALS
template<typename T>
complex<T> test0(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> tmp(-abs(last_z.real()), last_z.imag());

    const complex<T> new_z = tmp * tmp + c;

    return new_z;
}
template<typename T>
complex<T> test1(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> tmp(last_z.real(), -abs(last_z.imag()));

    const complex<T> new_z = tmp * tmp + c;

    return new_z;
}
template<typename T>
complex<T> test2(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    //For this the bailout can be 2 and the name of this fractal will be "Lakes" because the region in the top right corner reminds me of lakes
    const complex<T> tmp(last_z.real() * last_z.real() - last_z.imag() * last_z.imag(), T{2} * last_z.real() * last_z.real());

    const complex<T> new_z = tmp + c;

    return new_z;
}
template<typename T>
complex<T> test3(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test4(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test5(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test6(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test7(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test8(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
template<typename T>
complex<T> test9(const complex<T>& last_z, const complex<T>& c, [[maybe_unused]] const vector<complex<T>>& history, [[maybe_unused]] const vector<complex<T>>& extra_params, [[maybe_unused]] const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}