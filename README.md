cosmocalc
========
cosmocalc does basic cosmological calculations.

Matthew R. Becker, UChicago & Stanford/SLAC
Copyright (c) 2013-2014, New BSD License

USE THIS CODE AT YOUR OWN RISK! 

I have done some debugging, but I cannot promise 
that everything is correct. If you find any bugs, 
please let me know!

Matthew R. Becker  
becker.mr [at] gmail [dot] com 

Installation
--------------

You will need to install the GNU Scientific Library. Once you do 
that follow the instructions below.

C coders: You can compile static or shared libs. Also, if you look in examples/ you can see how to directly link the code.

Python Peeps: The entire C code has been wrapped by SWIG. To install, you 
should just be able to type
   
    python setup.py install

This assumes of course that your GSL is in a reasonable spot. If not, you 
will have to modify setup.py to include the paths to the headers and compiled 
libraries for GSL. (You can see this with the `/opt/local/lib options` in the 
script already. Just add your paths here. Note that environment variables are 
not resolved by python for whatever reason, so things like `$HOME/lib` will have 
to be resolved by hand. If you know how to do this, let me know!) 

Usage
---------

Units: I use h-inverse units, with Mpc/h and Msun/h where needed. Also, 
I work with scale factor as opposed to redshift, so be careful.

C coders: See the example directory. Note that cosmoData is a global struct with 
the cosmological parameters. When you want to change them, you need to increment 
field cosmoData.cosmoNum.

Python Peeps: It works like this

```python
import cosmocalc  
cd = {"om":0.3,"ob":0.045,"ol":0.7, \
      "ok":0.0,"h" :0.7,"s8":0.8,"ns":0.95, \
      "w0":-1.0,"wa":0.0}
cosmocalc.set_cosmology(cd)
print cosmocalc.comvdist(0.5)
```

The above code will print the comoving distance to scale factor = 0.5. Note that the 
cosmology is a *global* variable in C. Thus you can only have one cosmology working 
at a time. To set a new cosmology, simply change the dictionary and then call 
cosmocalc.set_cosmology().

The documentation is complete in the sense that in ipython, simply do this 

```python
>>> cosmocalc.comvdist?
```
to see what each function does.


How It Works
------------------

So the C code internally typically computes a given quantity as a function of 
scale factor, k, or whatever, and then builds a spline. The splines are built 
on-the-fly as needed. Thus the first call to a given quantity is slower than 
subsequent calls.

There is a single global variable (cosmoData) for the cosmology and the splines 
are rebuilt when the field cosmoData.cosmoNum is changed. The python interface 
does this automatically, but it has to be done by hand in the C code. 

The integrations are all done with various adaptive integrators in GSL. The 
details are not too important, except to say that the integrations are good to 
at least few percent (and a lot are much, much better). 

New BSD License
------------------------

Copyright (c) 2013-2014, Matthew R. Becker  
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
