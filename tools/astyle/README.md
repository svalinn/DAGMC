google-astyle
=============

Cpp source code formatter for Google C++ Style

[Astyle](http://astyle.sourceforge.net/astyle.html) is a free, fast and small
automatic formatter for C, C++, C+, and Java source code.

This fork, however, only aims at formatting C++ code according to [Google C++
Style Guide](http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml).

google-astyle code base is google-astyle clean, meaning you should always run
google-astyle on top of google-astyle source code before committing.

The command line flags we are using are:

```
astyle_google --style=linux               \
              --indent=spaces=2           \
              --pad-oper                  \
              --pad-header                \
              --unpad-paren               \
              --convert-tabs              \
              --align-pointer=type        \
              --suffix=none               \
              --recursive                 \
              --formatted                 \
              --lineend=linux             \
              --exclude=gtest             \
              *.cc *.cpp *.h *.hh *.hpp
```
