TO FINISH






/*
Chapter 1: what is a fractal? what does this program do?

A fractal, in mathematics, is a geometric shape containing detailed structure at arbitrarily small scales. (from Wikipedia)
In more informal terms, fractals are shapes with have details even if you zoom in on a portion of it 10000x times. And they have details at any
level of zoom. Their contours are usually very "rought", not straight, because the more you zoom, the more roughness you'll find.
Kind of.

The fractals that this program draws come from iterating complex functions of complex variables.
You have a function, f(z), and you repeatedly iterate over it like this
z = z0;
for(i = 0; i < 1000; ++i)
    z = f(z);
This simple loop sets z to f(z0), then f(f(z0)), then f(f(f(z0))) and so on.
z0 is said to be the initial condition of the iteration.

Usually the functions of the fractal depend on many parameters but the most common ones are:
- z: the value of f at the previous iteration
- c: the initial condition
For example, the famous Mandelbrot Set is defined as
f(z, c) = z^2 + c
and you compute the iterations of this function like this
z = 0;
c = (point in the complex plane);
for(i = 0; i < 1000; ++i)
    z = z^2 + c;

In order to draw the pretty and colorful images that you see around the internet (if you don't know what I'm talking about, go to youtube and
search for "fractal zoom"), you use a loop similar to this one, where the initial condition is set by sampling many different points in the complex
plane:
for(imag = -2; imag <= 2; imag += 0.01)
    for(real = -2; real <= 2; real += 0.01)
        iter_count = 0;
        z = 0;
        c = (real, imag);

        while(abs(z) < 2)
            z = z^2 + c;
            ++iter_count;
        
        image[x][y] = color(iter_count);
This simple loop iterates over a grid of points in the complex plane, iterates the Mandelbrot function until a condition is met and then colors
a pixel of the image with a color that depends on the amount of iterations it took to "diverge".

Note: for some points, the condition in the while loop is always true because those points never diverse, therefore it's usually replaced with
        [...]
        while(abs(z) < 2 && iter_count < max_iter)
            [...]
*/

/*
The problem that I'm facing arises from the fact that I want to have lots of options to tweak lots of core functions of the rendering process, the
most critical part of the program.
Lots of parameters are not known at compile time, but after the rendering is launched, those parameters stay constant.

In the beginning I had lots of duplicate code, I had a rendering function for each of the fractals. They were all very similar, the only difference
was the formula with which the next iteration was computed. Adding a new fractal meant lots of copy-pasting.
I then decided to try function pointers and write a single rendering loop, with the corresponding function for the next iteration called through
the pointer. This way I only had to set at runtime the pointer, depending on what the user chose to render and then launch the rendering process.
The result was a program ~10x slower.
I then tried using templates, with which I pass a function pointer as a template parameter and with this method the program is about as fast as
writing a separate rendering function for each fractal.

Wanting to tweak lots of core functions, though, means that I have more than 1 degree of freedom. I would like to maintain the performance of my
program as it is right now (and maybe improve it), but I can't find a decent and maintainable way to do it.
Passing the core parameters as a template means 

Below you can find more details of my implementation
*/

/*
I (kind-of) have a while loop, that I use to render the image, structured in the following way:
*/

template<fractal_fn_ptr_t fractal_fn>
void renderer(...){
    //[...]

    complex z, c;

    for(x = startX; x < endX; ++x){
        for(y = startY; y < endY; ++y){
            //Initialize variables z and c
            init_zc(z, c, x, y, ...);

            //Main loop
            iter_count = 0;
            while(iter_count < max_iter && (stop_iterating == false)){
                //Perform a single iteration dictacted by the chosen fractal function
                z = (*fractal_fn)(z, c, ...);

                if(exit_condition())
                    stop_iterating = true;
                
                ++iter_count;
            }

            image[x][y] = compute_color(iter_count, ...);
        }
    }
}

/*
The problem is that if I pass "fractal_fn" as a parameter to the function "renderer" the entire program becomes 10x slower.
As of now, since the user has a flag to choose which fractal to render, there is a big switch case in a function that does something like this:
*/
renderer_fn_ptr_t get_renderer_ptr(fractal_type){
    auto ret_ptr;

    switch(fractal_type){
        default:
        case fractal1: ret_ptr = &renderer<fractal1_fn>; break;
        case fractal2: ret_ptr = &renderer<fractal2_fn>; break;
        case fractal3: ret_ptr = &renderer<fractal3_fn>; break;
        case fractal4: ret_ptr = &renderer<fractal4_fn>; break;
        //[...]
    }

    return ret_ptr;
}
/*
A bunch of threads are then spawned, all of which call the same function through the pointer returned by "get_renderer_ptr".
*/

/*
Also, the main loop is written like this, therefore it's known at compile time as of now.
I already have another type of inner loop that is faster than this one in certain situations.
I want to add a flag that allows the user to choose what kind of inner loop to use without losing performance.
*/

/*
Also, as of now, the "exit_condition()" is just a simple comparison, known at compile time, but in the future I want to be able to
add a flag to the program that allows changing what kind of comparison is performed to exit.
*/