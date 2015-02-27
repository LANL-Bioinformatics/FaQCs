#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )
#exec >  >(tee install.log)
#exec 2>&1

cd $rootdir

install_perl_parallel_forkmanager()
{
echo "------------------------------------------------------------------------------
               Installing Perl Module Parallel-ForkManager-1.03
------------------------------------------------------------------------------
"
tar xvzf Parallel-ForkManager-1.03.tar.gz
cd Parallel-ForkManager-1.03
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/
cd $rootdir
rm -rf Parallel-ForkManager-1.03
echo "
------------------------------------------------------------------------------
                        Parallel-ForkManager-1.03 Installed
------------------------------------------------------------------------------
"
}

install_perl_string_approx()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module String-Approx-3.27
------------------------------------------------------------------------------
"
tar xvzf String-Approx-3.27.tar.gz
cd String-Approx-3.27
perl Makefile.PL 
make
cp -fR blib/lib/* $rootdir/
mkdir -p $rootdir/auto
mkdir -p $rootdir/auto/String
mkdir -p $rootdir/auto/String/Approx
cp -fR blib/arch/auto/String/Approx/Approx.* $rootdir/auto/String/Approx/
cd $rootdir
rm -rf String-Approx-3.27
echo "
------------------------------------------------------------------------------
                        String-Approx-3.27 Installed
------------------------------------------------------------------------------
"
}

checkPerlModule()
{
   perl -e "use lib \"$rootdir/lib\"; use $1;"
   return $?
}

checkSystemInstallation()
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}


if ( checkPerlModule Parallel::ForkManager )
then
  echo "Perl Parallel::ForkManager is found"
else
  echo "Perl Parallel::ForkManager is not found"
  install_perl_parallel_forkmanager
fi

if ( checkPerlModule String::Approx )
then
  echo "Perl String::Approx is found"
else
  echo "Perl String::Approx is not found"
  install_perl_string_approx
fi

if ( checkSystemInstallation R )
then
{
    echo "R is found."
}
else
{
    echo "R is not found"
    echo "Please install R from http://cran.r-project.org/";
    exit 1
}
fi

echo "

All done!

"
