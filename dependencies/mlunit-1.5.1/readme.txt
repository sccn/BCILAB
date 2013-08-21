mlUnit - Testing Framework for Matlab .m language.

README

===========================================================================

ABOUT

mlUnit is a testing framework for the Matlab .m language, considering the 
patterns of the xUnit family.

This software and all associated files are released unter the GNU General 
Public License (GPL) as published by the Free Software Foundation (see 
LICENSE file).

===========================================================================

INSTALLATION

1. Unzip mlunit.zip to $HOME.
2. Change to directory in MATLAB:

   >> cd $HOME/mlunit/src

3. Run install method:

   >> install(mlunit);

4. Add directory to MATLAB path:

   >> addpath('$HOME/mlunit/src');

===========================================================================

HOW TO TEST

As an example a test for the built-in sin function is written:

1. Create a new directory @test_sin: 

   >> mkdir @test_sin
   >> cd @test_sin

2. Create a new .m file test_sin.m (the constructor):

   >> edit test_sin.m

3. Add the following lines to test_sin.m:

   function self = test_sin(name)

   tc = test_case(name);
   self = class(struct([]), 'test_sin', tc);

4. Create a new file test_null.m (the first test) and add the following 
   lines:

   function self = test_null(self)

   assert_equals(0, sin(0));

5. Run the test:

   >> cd ..
   >> runner = text_test_runner(1, 1);
   >> loader = test_loader;
   >> run(runner, load_tests_from_test_case(loader, 'test_sin'));

   You should see something like this:

   .
   ----------------------------------------------------------------------
   Ran 1 test in 0.010s
   
   OK

6. Add more tests, e.g. test_sin_cos.m:

   function self = test_sin_cos(self)

   assert_equals(cos(0), sin(pi/2));

7. Run the tests:

   >> run(runner, load_tests_from_test_case(loader, 'test_sin'));
  
   You should see something like this:

   ..
   ----------------------------------------------------------------------
   Ran 2 tests in 0.010s
   
   OK

8. Try other parameters of text_test_runner:

   >> runner = text_test_runner(1, 2);
   >> run(runner, load_tests_from_test_case(loader, 'test_sin'));

   You should see something like this:

   test_null(test_sin) ... OK
   test_sin_cos(test_sin) ... OK

   ----------------------------------------------------------------------
   Ran 2 tests in 0.020s

   OK

9. Read more in the help texts of the classes, e.g.

   >> help text_test_runner

===========================================================================

MLUNIT TESTS

As mlUnit was developed test-driven, there are a number of tests in the 
test directory, which can be run by

>> addpath('$MLUNIT\mlunit\test');
>> mlunit_test;

===========================================================================

MORE DOCUMENTATION

More documentation can be found online at <http://mlunit.sourceforge.net>.

===========================================================================

QUESTIONS, COMMENTS, BUGS

If you have a question, a comment or a bug report, please send an email to 
thomas@dohmke.de. Usually I try to answer within 24 hours. Spam and other 
crap goes automatically to the trash bin, so please use a precise subject.