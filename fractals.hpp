#ifndef FRACTALS_HPP_INCLUDED
#define FRACTALS_HPP_INCLUDED

#include <png++/png.hpp>
#include <vector>
#include <cmath>
#include <complex>

#include <structs.h>
#include "misc_functions.hpp"

using namespace std;

extern atomic<int> runningThreads;

inline long double real_from_posX(const size_t posX, const size_t width, const size_t squareScale, const fractalsettings_t fset){
    return 4.0l * ((long double)posX - (long double)width / 2.0l)  / ((long double)squareScale * fset.scalingFactor) + fset.offset_re;
}

inline long double imag_from_posY(const size_t posY, const size_t height, const size_t squareScale, const fractalsettings_t fset){
    return 4.0l * ((long double)posY - (long double)height / 2.0l) / ((long double)squareScale * fset.scalingFactor) + fset.offset_im;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MANDELBROT SET
//z(n+1) = z(n)^2 +c
void mandelbrot(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    //z = re + i*im
    //c = cre + i*cim
    //z^2 + c = re^2 - im^2 + cre + i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = 2*re*im + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = real_from_posX(x, width, squareScale, fsettings);
            im = imag_from_posY(y, height, squareScale, fsettings);
            re_2 = re * re;
            im_2 = im * im;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < fsettings.max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate z^2 + c with the values calculated in the previous iteration
                im = 2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, re, im, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///TIPPETS MANDELBROT SET
//z(n+1) = z(n)^2 +c but wrongly implemented
/*
Normal:
im = 2 * re * im + cim;
re = re_2 - im_2 + cre;

Tippets:
re = re_2 - im_2 + cre;
im = 2 * re * im + cim;
*/
void tippets_mandelbrot(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    //z = re + i*im
    //c = cre + i*cim
    //z^2 + c = re^2 - im^2 + cre + i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = 2*re*im + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = real_from_posX(x, width, squareScale, fsettings);
            im = imag_from_posY(y, height, squareScale, fsettings);
            re_2 = re * re;
            im_2 = im * im;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < fsettings.max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate tippets mandelbrot set
                re = re_2 - im_2 + cre;
                im = 2 * re * im + cim;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, re, im, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///BURNING SHIP
//z(n+1) = (abs(re(z(n)) + i*abs(im(z(n)))) ^ 2 + c
void burningShip(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    //z = re + i*im
    //c = cre + i*cim
    //(abs(re) + i*abs(im))^2 + c = re^2 - im^2 + cre + abs(i*2*re*im) + i*cim;
    //re = re^2 - im^2 + cre
    //im = abs(2*re*im) + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = real_from_posX(x, width, squareScale, fsettings);
            im = imag_from_posY(y, height, squareScale, fsettings);
            re_2 = re * re;
            im_2 = im * im;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < fsettings.max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate burning ship with the values calculated in the previous iteration
                im = 2 * abs(re * im) + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, re, im, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MANDELBAR SET
//z(n+1) = conj(z(n))^2 +c
void mandelbar(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    //z = re + i*im
    //c = cre + i*cim
    //conj(z)^2 + c = re^2 - im^2 + cre - i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = -2*re*im + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = real_from_posX(x, width, squareScale, fsettings);
            im = imag_from_posY(y, height, squareScale, fsettings);
            re_2 = re * re;
            im_2 = im * im;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < fsettings.max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate conj(z)^2 + c with the values calculated in the previous iteration
                im = -2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, re, im, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///SPADE FRACTAL
//z(n+1) = z(n)^z(n) + z(n)/c
void spadefract(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //complex number a + ib = rho * e^(i*phi) (Euler formula)
    //rho = sqrt(a^2 + b^2), phi = atan2(b/a)

    //z(n+i) = z(n)^z(n) + z(n)/c
    //z  = zre  + i*zim  = zrho  * e^(i*zphi)
    //c  = cre  + i*cim
    //zz = zzre + i*zzim = zzrho * e^(i*zzphi)

    /*The following part demonstrates how to calculate z^z
    z^z = e^ln(z^z) = e^(z*ln(z))

    ln(z) = v = vre + i*vim
    ln(z) = ln(zrho * e^(i*zphi)) = ln(zrho) + ln(e^(i*zphi)) = ln(zrho) + i*zphi*ln(e) = ln(zrho) + i*zphi
    vre = ln(zrho) = ln(sqrt(zre^2 + zim^2))
    vim = zphi = atan2(zim/zre)

    z*ln(z) = z * v = (zre + i*zim) * (vre + i*vim) = w
    wre = zre * vre - zim * vim
    wim = zre * vim + zim * vre
    zzrho * e^(i*zzphi) = e^(wre + i*wim)
    e^(ln(zzrho) + i*zzphi) = e^(wre + i*wim) -> ln(zzrho) + i*zzphi = wre + i*wim
    ln(zzrho) = wre -> zzrho = e^wre
    zzphi = wim

    zzre = zzrho * cos(tan(zzphi))
    zzim = zzrho * sin(tan(zzphi))

    The following part demonstrates how to calculate z/c
    zc = zcre + i*zcim = z/c
    (zre + i*zim) / (cre + i*cim) =
    = ((zre + i*zim) * (cre - i*cim)) / ((cre + i*cim) * (cre - i*cim)) =
    = (zre * cre - zim * cim + i*(zre * cim + zim * cre)) / (cre^2 - cim^2)
    zcre = (zre * cre - zim * cim) / (cre^2 - cim^2)
    zcim = (zre * cim + zim * cre) / (cre^2 - cim^2)
    */

    //These are all the variables which represent all the complex numbers required for this fractal
    //The first two of each group are the cartesian representation, the other two are the polar representation
    long double zre   = 0.0l;
    long double zim   = 0.0l;
    long double zrho  = 0.0l;
    long double zphi  = 0.0l;

    long double zzre  = 0.0l;
    long double zzim  = 0.0l;
    long double zzrho = 0.0l;
    long double zzphi = 0.0l;

    long double cre   = 0.0l;
    long double cim   = 0.0l;

    long double vre   = 0.0l;
    long double vim   = 0.0l;

    long double wre   = 0.0l;
    long double wim   = 0.0l;

    long double zcre  = 0.0l;
    long double zcim  = 0.0l;

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            zre = real_from_posX(x, width, squareScale, fsettings);
            zim = imag_from_posY(y, height, squareScale, fsettings);
            zrho = sqrt(zre * zre + zim * zim);
            zphi = atan2(zim, zre);
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = zre;
                cim = zim;
            }

            int i = 0;
            while((i < fsettings.max_iter) && zrho < fsettings.bailout) {
                //Calculate z^z, all the calculations are explained earlier
                vre = log(zrho);
                vim = zphi;
                wre = zre * vre - zim * vim;
                wim = zre * vim + zim * vre;
                zzrho = exp(wre);
                zzphi = wim;
                zzre = zzrho * cos(tan(zzphi));
                zzim = zzrho * sin(tan(zzphi));
                //Calculate z/c, all the calculations are explained earlier
                zcre = (zre * cre - zim * cim) / (cre * cre - cim * cim);
                zcim = (zre * cim + zim * cre) / (cre * cre - cim * cim);

                //Assign to z the result of (z^z) + (z/c)
                zre = zzre + zcre;
                zim = zzim + zcim;

                //Calculate zrho and zphi for the check and the next iteration
                zrho = sqrt(zre * zre + zim * zim);
                zphi = atan2(zim, zre);
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, zre, zim, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MAGNET FRACTAL TYPE I
//z(n+1) = ((z(n)^2 + c - 1) / (2*z(n) + c - 2))^2
void magnet_type1(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    complex<long double> z{};
    complex<long double> c{};
    complex<long double> numerator{};
    complex<long double> denominator{};

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.juliaMode){
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

            int i = -1;
            while((i < fsettings.max_iter) && z.real() * z.real() + z.imag() * z.imag() < fsettings.bailout) {
                //z(n+1) = ((z(n)^2 + c - 1) / (2*z(n) + c - 2))^2
                numerator = z*z + c - complex<long double>{1, 0};
                denominator = complex<long double>{2, 0} * z + c - complex<long double>{2, 0};
                z = (numerator * numerator) / (denominator * denominator);
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MAGNET FRACTAL TYPE II
//z(n+1) = ((z(n)^3 + 3*(c - 1) * z(n) + (c - 1)*(c - 2)) / (3z(n)^2 + 3*(c - 2)*z(n) + (c - 1)*(c - 2) + 1))^2
void magnet_type2(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    complex<long double> z{};
    complex<long double> c{};
    complex<long double> numerator{};
    complex<long double> denominator{};

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.juliaMode){
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

            int i = -1;
            while((i < fsettings.max_iter) && z.real() * z.real() + z.imag() * z.imag() < fsettings.bailout) {
                //z(n+1) = ((z(n)^3 + 3*(c - 1) * z(n) + (c - 1)*(c - 2)) / (3z(n)^2 + 3*(c - 2)*z(n) + (c - 1)*(c - 2) + 1))^2
                numerator = z*z*z + complex<long double>{3, 0} * (c - complex<long double>{1, 0}) * z + (c - complex<long double>{1, 0}) * (c - complex<long double>{2, 0});    //(z(n)^3 + 3*(c - 1)*z(n) + (c - 1)*(c - 2))
                denominator = complex<long double>{3, 0} * z*z + complex<long double>{3, 0} * (c - complex<long double>{2, 0}) * z + (c - complex<long double>{1, 0}) * (c - complex<long double>{2, 0}) + complex<long double>{1, 0};  //3z(n)^2 + 3*(c - 2)*z(n) + (c - 1)*(c - 2) + 1

                z = (numerator * numerator) / (denominator * denominator);
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///CACTUS FRACTAL
//z(n+1) = z(n)^3 + (c - 1) * z(n) - c
void cactus(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    complex<long double> z{};
    complex<long double> c{};

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.juliaMode){
                z.real(real_from_posX(x, width, squareScale, fsettings));
                z.imag(imag_from_posY(y, height, squareScale, fsettings));
                c.real(fsettings.julia_re);
                c.imag(fsettings.julia_im);
            } else {
                /*
                By initializing z to zero we get a similar but different fractal
                z.real(0);
                z.imag(0);
                */

                z.real(real_from_posX(x, width, squareScale, fsettings));
                z.imag(imag_from_posY(y, height, squareScale, fsettings));
                c.real(real_from_posX(x, width, squareScale, fsettings));
                c.imag(imag_from_posY(y, height, squareScale, fsettings));
            }

            int i = -1;
            while((i < fsettings.max_iter) && z.real() * z.real() + z.imag() * z.imag() < fsettings.bailout) {
                //z(n+1) = z(n)^3 + (c - 1) * z(n) - c
                z = z*z*z + (c - complex<long double>{1, 0}) * z - c;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///ZUBIETA FRACTAL
//I don't know the real name of this fractal, I stole it from here http://paulbourke.net/fractals/Zubieta/
//z(n+1) = z(n)^2 + c / z(n)
void zubieta(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    /*
    z(n+1) = z(n)^2 + c / z(n)

    z = zre + i*zim
    new_z = new_zre + i*new_zim
    c = cre + i*cim
    cz = czre + i*czim                (stores the result of c/z)

    zre_2 = zre*zre   <-- updated at the end of the loop of the previous iteration
    zim_2 = zim*zim   <-- updated at the end of the loop of the previous iteration

    z^2 = zre^2 - zim^2 + i*2*zre*zim

    c/z = (cre + i*cim)/(zre + i*zim) =   (multiply by z* / z*)
    = (cre + i*cim)*(zre - i*zim)/(zre^2 + zim^2) =   (develop product)
    = (cre*zre - i*cre*zim + i*cim*zre + cim*zim)/(zre^2 + zim^2)
    czre = (cre*zre + cim*zim) / (zre^2 + zim^2)
    czim = (-i*cre*zim + i*cim*zre) / (zre^2 + zim^2)

    new_z = z^2 + z/c
    new_zre = zre^2 - zim^2 + czre = zre^2 - zim^2 + (cre*zre + cim*zim) / (zre^2 + zim^2)
    new_zim = 2*zre*zim + czim = 2*zre*zim + (-cre*zim + cim*zre) / (zre^2 + zim^2)
    */

    long double zre     = 0.0l;
    long double zim     = 0.0l;
    long double new_zre = 0.0l;
    long double new_zim = 0.0l;
    long double zre_2   = 0.0l;    //zre^2
    long double zim_2   = 0.0l;    //zim^2
    long double cre     = 0.0l;
    long double cim     = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            zre = real_from_posX(x, width, squareScale, fsettings);
            zim = imag_from_posY(y, height, squareScale, fsettings);
            zre_2 = zre * zre;
            zim_2 = zim * zim;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = zre;
                cim = zim;
            }

            int i = 0;
            while((i < fsettings.max_iter) && (zre_2 + zim_2) < bailout_2) {
                //Calculate z^2 + c/z
                new_zre = zre_2 - zim_2 + (cre*zre + cim*zim) / (zre_2 + zim_2);
                new_zim = 2*zre*zim + (-1*cre*zim + cim*zre) / (zre_2 + zim_2);

                //Update z
                zre = new_zre;
                zim = new_zim;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                zre_2 = zre * zre;
                zim_2 = zim * zim;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, zre, zim, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///ZUBITHETA FRACTAL
//Because after eta comes theta. This is very similar to the Zubieta fractal
//z(n+1) = z(n)^2 + z(n) / c
void zubitheta(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    /*
    z(n+1) = z(n)^2 + z(n) / c

    z = zre + i*zim
    new_z = new_zre + i*new_zim
    c = cre + i*cim
    zc = zcre + i*zcim                (stores the result of z/c)

    zre_2 = zre*zre   <-- updated at the end of the loop of the previous iteration
    zim_2 = zim*zim   <-- updated at the end of the loop of the previous iteration

    z^2 = zre^2 - zim^2 + i*2*zre*zim

    z/c = (zre + i*zim)/(cre + i*cim) =   (multiply by c* / c*)
    = (zre + i*zim)*(cre - i*cim)/(cre^2 + cim^2) =   (develop product)
    = (zre*cre + zim*cim + i*(zre*(-cim) + zim*cre)) / (cre^2 + cim^2)
    zcre = (zre*cre + zim*cim) / (cre^2 + cim^2)
    zcim = (zre*(-cim) + zim*cre)) / (cre^2 + cim^2)

    new_z = z^2 + c/z
    new_zre = zre^2 - zim^2 + zcre = zre^2 - zim^2 + (zre*cre + zim*cim) / (cre^2 + cim^2)
    new_zim = 2*zre*zim + zcim = 2*zre*zim + (zre*(-cim) + zim*cre)) / (cre^2 + cim^2)
    cre^2 + cim^2 = cnorm
    */

    long double zre     = 0.0l;
    long double zim     = 0.0l;
    long double new_zre = 0.0l;
    long double new_zim = 0.0l;
    long double zre_2   = 0.0l;    //zre^2
    long double zim_2   = 0.0l;    //zim^2
    long double cre     = 0.0l;
    long double cim     = 0.0l;
    long double cnorm   = 0.0l;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            zre = real_from_posX(x, width, squareScale, fsettings);
            zim = imag_from_posY(y, height, squareScale, fsettings);
            zre_2 = zre * zre;
            zim_2 = zim * zim;
            if(fsettings.juliaMode){
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                cre = zre;
                cim = zim;
            }
            cnorm = cre * cre + cim * cim;

            int i = 0;
            while((i < fsettings.max_iter) && (zre_2 + zim_2) < bailout_2) {
                //Calculate z^2 + c/z
                new_zre = zre_2 - zim_2 + (zre*cre + zim*cim) / cnorm;
                new_zim = 2*zre*zim + (zre*(-1)*cim + zim*cre) / cnorm;

                //Update z
                zre = new_zre;
                zim = new_zim;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                zre_2 = zre * zre;
                zim_2 = zim * zim;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, zre, zim, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///LOGISTIC MAP FRACTAL
//z(n+1) = c * z(n) * (1-z(n))
void logistic_map(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = fsettings.bailout * fsettings.bailout;

    /*
    z(n+1) = c * z(n) * (1-z(n))

    z = zre + i*zim
    new_z = new_zre + i*new_zim
    c = cre + i*cim

    (cre + i*cim) * (zre + i*zim) * (1 - zre - i*zim)

    Develop the last product first and call this result w
    w = (zre + i*zim)*((1-zre) + i*(-zim)) =
    = zre*(1-zre) + zim*zim + i*(zre*(-zim) + zim*(1-zre)) =
    = zre - zre*zre + zim*zim + i*(-1*zre*zim + zim - zim*zre)
    wre = zre - zre*zre + zim*zim
    wim = -1*zre*zim + zim - zim*zre

    Good, now time to compute c*w
    new_z = (cre + i*cim) * (wre + i*wim) =
    = cre*wre - cim*wim + i*(cre*wim + cim*wre)
    */

    long double zre     = 0.0l;
    long double zim     = 0.0l;
    long double new_zre = 0.0l;
    long double new_zim = 0.0l;
    long double zre_2   = 0.0l;     //zre^2
    long double zim_2   = 0.0l;     //zim^2
    long double cre     = 0.0l;
    long double cim     = 0.0l;
    long double wre     = 0.0l;     //real(z*(1-z))
    long double wim     = 0.0l;     //imag(z*(1-z))

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.juliaMode){
                zre = real_from_posX(x, width, squareScale, fsettings);
                zim = imag_from_posY(y, height, squareScale, fsettings);
                cre = fsettings.julia_re;
                cim = fsettings.julia_im;
            } else {
                zre = 0.5l;
                zim = 0;
                cre = real_from_posX(x, width, squareScale, fsettings);
                cim = imag_from_posY(y, height, squareScale, fsettings);
            }
            zre_2 = zre * zre;
            zim_2 = zim * zim;

            int i = -1;
            while((i < fsettings.max_iter) && (zre_2 + zim_2) < bailout_2) {
                //Calculate w = z*(1-z)
                wre = zre - zre_2 + zim_2;
                wim = -1*zre*zim + zim - zim*zre;

                //Calculate new_z = c*w
                new_zre = cre*wre - cim*wim;
                new_zim = cre*wim + cim*wre;

                //Update z
                zre = new_zre;
                zim = new_zim;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                zre_2 = zre * zre;
                zim_2 = zim * zim;
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, zre, zim, fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///BINOCULAR FRACTAL (m = 4)
//z(n+1) = (z(n)^m + z(n)^2 + 1) / (2*z(n)^(m-1) - c + 1)
void binocular_m4(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    complex<long double> z{};
    complex<long double> c{};

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            if(fsettings.juliaMode){
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


            int i = -1;
            while((i < fsettings.max_iter) && z.real() * z.real() + z.imag() * z.imag() < fsettings.bailout) {
                //z(n+1) = z(n+1) = (z(n)^4 + z(n)^2 + 1) / (2*z(n)^3 - c + 1)
                z = (pow(z, 4) + pow(z, 2) + complex<long double>{1, 0}) / (complex<long double>{2, 0} * pow(z, 3) - c + complex<long double>{1, 0});
                i++;
            }

            image_to_write[y][x] = compute_color(i, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.colorMode, palette);
        }
    }

    return;
}

///TODO
//Find other fractals

#endif // FRACTALS_HPP_INCLUDED
