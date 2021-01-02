#ifndef FRACTALS_HPP_INCLUDED
#define FRACTALS_HPP_INCLUDED

#include <png++/png.hpp>
#include <vector>
#include <atomic>
#include <cmath>

#include "misc_functions.hpp"

using namespace std;

extern atomic<int> runningThreads;

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MANDELBROT SET FUNCTIONS
//z(n+i) = z(n)^2 +c
/*
vector<vector<size_t>> mandelbrot_iter(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const double bailout,
    const size_t width, const size_t height) {

    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

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

    vector<vector<size_t>> return_iter(endY - startY, vector<size_t>(endX - startX, 0));

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate z^2 + c with the values calculated in the previous iteration
                im = 2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            //Store calculated amounts of iteration in return_iter vector of vectors
            return_iter[y - startY][x - startX] = i;
        }
    }

    runningThreads -= 1;
    return return_iter;
}
*/

void mandelbrot_color(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

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
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate z^2 + c with the values calculated in the previous iteration
                im = 2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, max_iter, re, im, bailout, color_mode, palette);
        }
    }

    runningThreads -= 1;
    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///BURNING SHIP FUNCTIONS
//z(n+1) = (abs(re(z(n)) + i*abs(im(z(n)))) ^ 2 + c
/*
vector<vector<size_t>> burningShip_iter(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const double bailout,
    const size_t width, const size_t height) {

    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

    //z = re + i*im
    //c = cre + i*cim
    //(abs(re) + i*abs(im)) ^ 2 + c = re^2 - im^2 + 2*i*abs(re)*abs(im) + c =
    //= re^2 - im^2 + 2*abs(re*im) + c

    //re = re^2 - im^2 + cre
    //im = 2*re*im + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    vector<vector<size_t>> return_iter(endY - startY, vector<size_t>(endX - startX, 0));

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate burning ship with the values calculated in the previous iteration
                im = 2 * abs(re * im) + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            //Store calculated amounts of iteration in return_iter vector of vectors
            return_iter[y - startY][x - startX] = i;
        }
    }

    runningThreads -= 1;
    return return_iter;
}
*/

void burningShip_color(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

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
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate burning ship with the values calculated in the previous iteration
                im = 2 * abs(re * im) + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, max_iter, re, im, bailout, color_mode, palette);
        }
    }

    runningThreads -= 1;
    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///MANDELBAR SET FUNCTIONS
//z(n+i) = conj(z(n))^2 +c
/*
vector<vector<size_t>> mandelbar_iter(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const double bailout,
    const size_t width, const size_t height) {

    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

    //z = re + i*im
    //c = cre + i*cim
    //conj(z)^2 + c = re^2 - im^2 + cre - i*2*re*im + i*cim;
    //re = re^2 - im^2 + cre
    //im = 2*re*im + cim

    long double re   = 0.0l;
    long double im   = 0.0l;
    long double re_2 = 0.0l;    //re^2
    long double im_2 = 0.0l;    //im^2
    long double cre  = 0.0l;
    long double cim  = 0.0l;

    vector<vector<size_t>> return_iter(endY - startY, vector<size_t>(endX - startX, 0));

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate z^2 + c with the values calculated in the previous iteration
                im = -2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            //Store calculated amounts of iteration in return_iter vector of vectors
            return_iter[y - startY][x - startX] = i;
        }
    }

    runningThreads -= 1;
    return return_iter;
}
*/

void mandelbar_color(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto squareScale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout_2 = bailout * bailout;

    //z = re + i*im
    //c = cre + i*cim
    //conj(z)^2 + c = re^2 - im^2 + cre - i*2*re*im + i*cim;
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
            re = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            im = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            re_2 = re * re;
            im_2 = im * im;
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = re;
                cim = im;
            }

            int i = 0;
            while((i < max_iter) && (re_2 + im_2) < bailout_2) {
                //Calculate z^2 + c with the values calculated in the previous iteration
                im = -2 * re * im + cim;
                re = re_2 - im_2 + cre;

                //Calculate re^2 and im^2 for next iteration and easier norm checking;
                re_2 = re * re;
                im_2 = im * im;
                i++;
            }

            image_to_write[y][x] = compute_color(i, max_iter, re, im, bailout, color_mode, palette);
        }
    }

    runningThreads -= 1;
    return;
}

///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
///SPADE FRACTAL FUNCTIONS
//z(n+i) = z(n)^z(n) + z(n)/c
/*
vector<vector<size_t>> spadefract_iter(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const double bailout,
    const size_t width, const size_t height) {

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

    //The following part demonstrates how to calculate z^z
    //z^z = e^ln(z^z) = e^(z*ln(z))

    //ln(z) = v = vre + i*vim
    //ln(z) = ln(zrho * e^(i*zphi)) = ln(zrho) + ln(e^(i*zphi)) = ln(zrho) + i*zphi*ln(e) = ln(zrho) + i*zphi
    //vre = ln(zrho) = ln(sqrt(zre^2 + zim^2))
    //vim = zphi = atan2(zim/zre)

    //z*ln(z) = z * v = (zre + i*zim) * (vre + i*vim) = w
    //wre = zre * vre - zim * vim
    //wim = zre * vim + zim * vre
    //zzrho * e^(i*zzphi) = e^(wre + i*wim)
    //e^(ln(zzrho) + i*zzphi) = e^(wre + i*wim) -> ln(zzrho) + i*zzphi = wre + i*wim
    //ln(zzrho) = wre -> zzrho = e^wre
    //zzphi = wim

    //zzre = zzrho * cos(tan(zzphi))
    //zzim = zzrho * sin(tan(zzphi))

    //The following part demonstrates how to calculate z/c
    //zc = zcre + i*zcim = z/c
    //(zre + i*zim) / (cre + i*cim) =
    //= ((zre + i*zim) * (cre - i*cim)) / ((cre + i*cim) * (cre - i*cim)) =
    //= (zre * cre - zim * cim + i*(zre * cim + zim * cre)) / (cre^2 - cim^2)
    //zcre = (zre * cre - zim * cim) / (cre^2 - cim^2)
    //zcim = (zre * cim + zim * cre) / (cre^2 - cim^2)

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

    vector<vector<size_t>> return_iter(endY - startY, vector<size_t>(endX - startX, 0));

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            zre = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            zim = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            zrho = sqrt(zre * zre + zim * zim);
            zphi = atan2(zim, zre);
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = zre;
                cim = zim;
            }

            int i = 0;
            while((i < max_iter) && zrho < bailout) {
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

            //Store calculated amounts of iteration in return_iter vector of vectors
            return_iter[y - startY][x - startX] = i;
        }
    }

    runningThreads -= 1;
    return return_iter;
}
*/

void spadefract_color(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write) {

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

    //The following part demonstrates how to calculate z^z
    //z^z = e^ln(z^z) = e^(z*ln(z))

    //ln(z) = v = vre + i*vim
    //ln(z) = ln(zrho * e^(i*zphi)) = ln(zrho) + ln(e^(i*zphi)) = ln(zrho) + i*zphi*ln(e) = ln(zrho) + i*zphi
    //vre = ln(zrho) = ln(sqrt(zre^2 + zim^2))
    //vim = zphi = atan2(zim/zre)

    //z*ln(z) = z * v = (zre + i*zim) * (vre + i*vim) = w
    //wre = zre * vre - zim * vim
    //wim = zre * vim + zim * vre
    //zzrho * e^(i*zzphi) = e^(wre + i*wim)
    //e^(ln(zzrho) + i*zzphi) = e^(wre + i*wim) -> ln(zzrho) + i*zzphi = wre + i*wim
    //ln(zzrho) = wre -> zzrho = e^wre
    //zzphi = wim

    //zzre = zzrho * cos(tan(zzphi))
    //zzim = zzrho * sin(tan(zzphi))

    //The following part demonstrates how to calculate z/c
    //zc = zcre + i*zcim = z/c
    //(zre + i*zim) / (cre + i*cim) =
    //= ((zre + i*zim) * (cre - i*cim)) / ((cre + i*cim) * (cre - i*cim)) =
    //= (zre * cre - zim * cim + i*(zre * cim + zim * cre)) / (cre^2 - cim^2)
    //zcre = (zre * cre - zim * cim) / (cre^2 - cim^2)
    //zcim = (zre * cim + zim * cre) / (cre^2 - cim^2)

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

    vector<vector<size_t>> return_iter(endY - startY, vector<size_t>(endX - startX, 0));

    //Calculate the iterations of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize doubles
            zre = 4.0l * ((long double)x - (long double)width / 2.0l)  / ((long double)squareScale * scalingFactor) + offset_re;
            zim = 4.0l * ((long double)y - (long double)height / 2.0l) / ((long double)squareScale * scalingFactor) + offset_im;
            zrho = sqrt(zre * zre + zim * zim);
            zphi = atan2(zim, zre);
            if(juliaMode){
                cre = julia_re;
                cim = julia_im;
            } else {
                cre = zre;
                cim = zim;
            }

            int i = 0;
            while((i < max_iter) && zrho < bailout) {
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

            image_to_write[y][x] = compute_color(i, max_iter, zre, zim, bailout, color_mode, palette);
        }
    }

    runningThreads -= 1;
    return;
}

///TODO
//Magnet fractal
//Find other fractals

#endif // FRACTALS_HPP_INCLUDED
