#ifndef FRACTALS_HPP_INCLUDED
#define FRACTALS_HPP_INCLUDED

#include <vector>
#include <cmath>
#include <complex>

#include "structs.hpp"
#include "misc_functions.hpp"

# define M_PIl          3.141592653589793238462643383279502884L

using namespace std;

inline long double real_from_posX(const size_t posX, const size_t width, const size_t squareScale, const fractalsettings_t fset);
inline long double imag_from_posY(const size_t posY, const size_t height, const size_t squareScale, const fractalsettings_t fset);

void init_fractal(const size_t& x, const size_t& y, const size_t& width, const size_t& height, const size_t& square_scale, complex<long double>& z, complex<long double>& c, vector<complex<long double>>& history, vector<complex<long double>>& extra_params, const fractalsettings_t& fset);

//Mandelbrot set
//z(n+1) = z(n)^2 +c
template<typename T>
complex<T> mandelbrot(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//TIPPETS MANDELBROT SET
//z(n+1) = z(n)^2 +c but wrongly implemented
template<typename T>
complex<T> tippets_mandelbrot(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//BURNING SHIP
//z(n+1) = (abs(re(z(n)) + i*abs(im(z(n)))) ^ 2 + c
template<typename T>
complex<T> burning_ship(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//MANDELBAR SET
//z(n+1) = conj(z(n))^2 +c
template<typename T>
complex<T> mandelbar(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//MAGNET FRACTAL TYPE I
//z(n+1) = ((z(n)^2 + c - 1) / (2*z(n) + c - 2))^2
template<typename T>
complex<T> magnet_type1(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//MAGNET FRACTAL TYPE II
//z(n+1) = ((z(n)^3 + 3*(c - 1) * z(n) + (c - 1)*(c - 2)) / (3z(n)^2 + 3*(c - 2)*z(n) + (c - 1)*(c - 2) + 1))^2
template<typename T>
complex<T> magnet_type2(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//CACTUS FRACTAL
//CACTUS0 FRACTAL
//z(n+1) = z(n)^3 + (c - 1) * z(n) - c
template<typename T>
complex<T> cactus(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//ZUBIETA FRACTAL
//z(n+1) = z(n)^2 + c / z(n)
template<typename T>
complex<T> zubieta(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//ZUBITHETA FRACTAL
//z(n+1) = z(n)^2 + z(n) / c
template<typename T>
complex<T> zubitheta(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//LOGISTIC MAP FRACTAL
//z(n+1) = c * z(n) * (1-z(n))
template<typename T>
complex<T> logistic_map(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//UNPOLISHED SQUARE FRACTAL
//Re(n) = (Re(n-1) - Im(n-1) )*|Im(n-1)| + Re_c
//Im(n) = (Re(n-1) + Im(n-1) )*|Re(n-1)| + Im_c
template<typename T>
complex<T> unpol_square(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//MOTH FRACTAL
//z(n+1) = (3*z^3 - 2*z - 1)/(z * c + 1)
template<typename T>
complex<T> moth(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//CAVE FRACTAL
//z(n+1) = (z(n)^3 + c)/(-2*z(n) + 1)
template<typename T>
complex<T> cave(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//WANKEL FRACTAL
//z(n+1) = (z(n)^3 + 1)/c
template<typename T>
complex<T> wankel(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//SEA ANGEL FRACTAL
//z(n+1) = (z(n)^4 + 3*z(n)^2 + c) / (5*z(n)^2 - 3*z(n) + 2)
template<typename T>
complex<T> sea_angel(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//SMITH CHART FRACTAL
//z(n+1) = (1+z(n)) / (1-z(n)) + c
template<typename T>
complex<T> smith(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//SPADE FRACTAL
//z(n+1) = z(n)^z(n) + z(n)/c
template<typename T>
complex<T> spadefract(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

//TEST FRACTALS
template<typename T>
complex<T> test0(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test1(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test2(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test3(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test4(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test5(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test6(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test7(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test8(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);
template<typename T>
complex<T> test9(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings);

/*
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TESTING FRACTALS
//The formula in this fractal is subject to change from time to time to be able to quickly test new fractals
void test(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    const static size_t historySize = 128;        //To hold z(n-1) and z(n-2)
    vector<complex<T>> z{};
    complex<T> new_z{};
    complex<T> c{};

    //Calculate the iterations of the pixels for each posize_t in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.julia_mode){
                z = vector<complex<T>>(historySize, complex<T>(real_from_posX(x, width, squareScale, fsettings), imag_from_posY(y, height, squareScale, fsettings)));
                c.real(fsettings.julia_re);
                c.imag(fsettings.julia_im);
            } else {
                z = vector<complex<T>>(historySize, complex<T>(0, 0));
                c.real(real_from_posX(x, width, squareScale, fsettings));
                c.imag(imag_from_posY(y, height, squareScale, fsettings));
            }

            size_t i = -1;
            while((i < fsettings.max_iter) && (z.front().real() * z.front().real() + z.front().imag() * z.front().imag()) < fsettings.bailout) {
                //z(n) = z(n-1)*z(n-2) + c
                new_z = accumulate(z.begin(), z.end(), complex<T>(1,0), multiplies<complex<T>>()) + c;
                z.pop_back();
                z.insert(z.begin(), new_z);
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.front().real(), z.front().imag(), fsettings.bailout, csettings.cmode, palette);
        }
    }

    return;
}

//The formula in this fractal is subject to change from time to time to be able to quickly test new fractals
void test2(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    complex<T> z{};
    complex<T> c{};

    //Calculate the iterations of the pixels for each posize_t in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.julia_mode){
                z.real(real_from_posX(x, width, squareScale, fsettings));
                z.imag(imag_from_posY(y, height, squareScale, fsettings));
                c.real(fsettings.julia_re);
                c.imag(fsettings.julia_im);
            } else {
                z.real(0);
                z.imag(0);
                c.real(real_from_posX(x, width, squareScale, fsettings));
                c.imag(imag_from_posY(y, height, squareScale, fsettings));
            }

            size_t i = -1;
            while((i < fsettings.max_iter) && (z.real() * z.real() + z.imag() * z.imag()) < fsettings.bailout) {
                //z(n) = z(n-1)*z(n-2) + c
                z = complex<T>(0.25, 0) * (complex<T>(2, 0) + z+z+z+z+z+z+z - (complex<T>(2, 0) + z+z+z+z+z) * cos(z * complex<T>(M_PIl, 0)));
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.cmode, palette);
        }
    }

    return;
}
*/

//TODO
//Find other fractals
//BINOCULAR FRACTAL (m = 3, 4, 5)
//z(n+1) = (z(n)^m + z(n)^2 + 1) / (2*z(n)^(m-1) - c + 1)
//https://www.reddit.com/r/fractals/comments/r80ahw/polynomial_quotient_fractal_and_julia_set_z_az_bz/?utm_medium=android_app&utm_source=share
//https://www.reddit.com/r/fractals/comments/reo9d0/factorized_polynomial_quotient_fractal_zoom_and/?utm_medium=android_app&utm_source=share

/*
Strange magnet, seems like an intermediate mandelbrot set, like one where the power is 1.5 or sth
complex<T> magnet_type_mandel(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings) {
    const complex<T> numerator = last_z * last_z + c - complex<T>(1, 0);
    const complex<T> denominator = complex<T>(2, 0) * last_z + c - complex<T>(2, 0);
    
    const complex<T> new_z = (numerator * numerator) / (denominator);

    return new_z;
}


Another fractal that came up while rewriting magnet type 2
Edit 43 seconds later: i found the error, the denominator is not squared
The structure is a bit phallic lol but I can't call this the dick fractal
complex<T> magnet_type_error(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings) {
    const complex<T> numerator = last_z*last_z*last_z + complex<T>(3, 0) * (c - complex<T>(1, 0)) * last_z + (c - complex<T>(1, 0)) * (c - complex<T>(2, 0));
    const complex<T> denominator = complex<T>(3, 0) * last_z*last_z + complex<T>(3, 0) * (c - complex<T>(2, 0)) * last_z + (c - complex<T>(1, 0)) * (c - complex<T>(2, 0)) + complex<T>(1, 0);
    
    const complex<T> new_z = (numerator * numerator) / (denominator);

    return new_z;
}

Buffalo fractal
z = z^2 - z + c;

Perpendicular mandelbrot
z = [abs(re(z))-i*im(z)]^2 + c

Other perpendicular fractals
https://www.deviantart.com/kosmic-stardust/art/Mandelbrot-ABS-Variations-Complete-Set-of-Formulas-487039852
Cubic variations
https://www.deviantart.com/kosmic-stardust/art/Cubic-Mandelbrot-ABS-Variations-Incomplete-487039945

Fractals from fractal sound explorer from here
https://codeparade.itch.io/fractal-sound-explorer
Feather fractal z = z^3/(1+z^2) + c (note: the z^2 might actually be a dot product or a cross product, it had a 0 symbol)
SFX fractal z*(z dot z) - z * (c 0 c) again, idk what these symbols mean
Chirikov map: y_new = y + c_y * sin(x), x_new = x + c_x * y_new
Henon map: https://en.wikipedia.org/wiki/H%C3%A9non_map
Duffing map: https://en.wikipedia.org/wiki/Duffing_map
Ikeda map: https://en.wikipedia.org/wiki/Ikeda_map

Frog fractal
z = sin(z^2) + c^z - c
with z0=0
bailout 128
source: https://www.youtube.com/watch?v=8d4XRauHlj8

z = 5*sin(5 * z^5 + 5) + 5^c
bailout 256
source: https://www.youtube.com/watch?v=lheRzRHUjsI

Implement one of these or all of these mandelbrot inversion formulas
https://www.youtube.com/watch?v=MULce2qryes

formula for quasi davisbrot:
f(z)=(z.y^3-z.x^2+x,-2z.xz.y+y)
formula for quasi davisbar:
f(z)=(z.y^3-z.x^2+x,2z.xz.y+y)
formula for davisbrot:
f(z)=(z.x^2-z.y^3+x,2z.xz.y+y)
formula for davisbar:
f(z)=(z.x^2-z.y^3+x,-2z.xz.y+y)
formula for perp. davisbrot:
f(z)=(z.x^3-z.y^3+x,2z.xz.y+y)
formula for quasi perp. davisbrot:
f(z)=(z.x^3-z.y^3+x,-2z.xz.y+y)
formula for burning davisbrot:
f(z)=(z.x^2-z.y+x,-2z.xz.y+y)
formula for q. burning davisbrot:
f(z)=(z.y-z.x^2+x,-2z.xz.y+y)
second formula for ALL:
D(z)=sqrt(z.x^2+z.y^2)
*/

/*
Template for new fractals
template<typename T>
complex<T> template(const complex<T>& last_z, const complex<T>& c, const vector<complex<T>>& history, const vector<complex<T>>& extra_params, const fractalsettings_t& fsettings) {
    const complex<T> new_z = last_z * last_z + c;

    return new_z;
}
*/

#include "fractals.ipp"

#endif // FRACTALS_HPP_INCLUDED
