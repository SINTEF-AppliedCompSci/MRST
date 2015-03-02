#!/usr/bin/pyton
import sys
import os

def configure_dir(rootdir,rootdir_codes,opt_dir,cpp_dir,dolink):
    mydir = rootdir_build +'/' + opt_dir + '/' + cpp_dir
    if(os.path.exists(mydir)==False):
        os.mkdir(mydir)
    os.chdir(mydir)
    command =  rootdir_codes + '/' + cpp_dir + '/configure'
    #command = 'eval ' + rootdir_codes + '/' + cpp_dir + '/configure'
    if opt_dir=='db':
        command = command + '  ' + db_conf
        if(dolink==True):
            command = command + ' ' + link_options_db
    elif  opt_dir=='opt':
        command = command + ' ' + release_conf
        if(dolink==True):
            command = command + ' ' + link_options_release
    else:
        print 'no such opt ' + opt_dir
           
    print "Command : " + command
    print "Dir " + os.getcwd()
    a= os.system(command)
    if a!=0:
        print 'Command failed :' + command
        print 'In dir ' + mydir
        sys.exit(3)
    return

def configure_dirs(rootdir_build, rootdir_codes,comopt,dir_in,dolink):
    for opt_dir in comopt:
        if(os.path.exists(rootdir_build + '/' + opt_dir)==False):
            os.mkdir(rootdir_build + '/' + opt_dir)
        for cpp_dir in dir_in:
            a=os.path.exists(rootdir_build + '/' + opt_dir + '/' +cpp_dir + '/config.status')
            if(a==False):
                configure_dir(rootdir_build,rootdir_codes,opt_dir,cpp_dir,dolink)
    return

def make_dirs(rootdir_build,comopt,dir_in,inst):
    for opt_dir in comopt:
        for cpp_dir in dir_in:
            mydir = rootdir_build +'/' + opt_dir + '/' + cpp_dir
            if(os.path.exists(mydir)==False):
                os.mkdir(mydir)            
            os.chdir(mydir) 
            a= os.system('make -j 4')
            if(a==0):
                if(inst==True):
                    b = os.system('make install')
                    if(b>0):
                        print 'Instalation failed ' + mydir
                        sys.exit(3)
            else:
                print 'Compilation failed in ' + mydir
                sys.exit(3)
    return
        
#first build all tings which need to be installed
myrootdir=os.getcwd() + '/' + 'codes'
rootdir_build= myrootdir + '/builds'
rootdir_codes= myrootdir
rootdir_agmg = '/home/hnil/PROGRAM/AGMG_3.1.2'
#rootdir_agmg = myrootdir + '/AGMG'

# set op configurations
link_options_db='LDFLAGS="-L' + rootdir_codes + '/inst/db/lib" CPPFLAGS="-I' + rootdir_codes + '/inst/db/include"'
link_options_release='LDFLAGS="-L' + rootdir_codes + '/inst/opt/lib" CPPFLAGS="-I' + rootdir_codes + '/inst/opt/include"'

if(os.path.exists(rootdir_agmg + '/SRC')):
   with_agmg="--with-agmg=" + rootdir_agmg + "/SRC "
else:
   with_agmg=""
   
db_conf="CC=gcc CXX=g++ F77=gfortran-4.4 FC=gfortran-4.4 --enable-static " + with_agmg + " --disable-shared CFLAGS='-g -Wall -Wextra  -ansi -pedantic -Wformat-nonliteral -Wcast-align -Wpointer-arith -Wbad-function-cast -Wmissing-prototypes -Wstrict-prototypes -Wmissing-declarations -Winline -Wundef -Wnested-externs -Wcast-qual -Wshadow -Wwrite-strings -Wchar-subscripts -Wredundant-decls -Wno-unknown-pragmas' CXXFLAGS='-ggdb3 -O0 -Wall -Wextra  -ansi -pedantic  -Wformat-nonliteral -Wcast-align -Wpointer-arith -Wmissing-declarations -Wundef -Wcast-qual -Wshadow -Wwrite-strings -Wchar-subscripts -Wredundant-decls -Wno-unknown-pragmas -D_GLIBCXX_DEBUG=1 -D_GLIBCXX_DEBUG_PEDANTIC=1' --prefix=" + rootdir_codes + "/inst/db/ INSTALL=\'install -C'"
release_conf="CC=gcc CXX=g++ F77=gfortran-4.4 FC=gfortran-4.4 --enable-static " + with_agmg + " --disable-shared CFLAGS='-DNDEBUG -O3 -Wall -Wextra  -ansi -pedantic -Wformat-nonliteral -Wcast-align -Wpointer-arith -Wbad-function-cast -Wmissing-prototypes -Wstrict-prototypes -Wmissing-declarations -Winline -Wundef -Wnested-externs -Wcast-qual -Wshadow -Wwrite-strings -Wchar-subscripts -Wredundant-decls -Wno-unknown-pragmas' CXXFLAGS='-DNDEBUG -O3 -Wall -Wextra  -ansi -pedantic  -Wformat-nonliteral -Wcast-align -Wpointer-arith -Wmissing-declarations -Wundef -Wcast-qual -Wshadow -Wwrite-strings -Wchar-subscripts -Wredundant-decls -Wno-unknown-pragmas' --prefix=" + rootdir_codes + "/inst/opt/ INSTALL='install -C'"
# define directoies for db and release 
#comopt=['db','opt']
comopt=['db','opt']
for optdir in comopt:
    mydir=rootdir_build + '/' + 'inst' + '/' + optdir
    a=os.path.exists(mydir)
    if(a==False):
        os.makedirs(mydir) 


if(os.path.exists(rootdir_build)==False):
    os.mkdir(rootdir_build)
#os.mkdir(rootdir_condes)
# define directories which whould be installed
if(os.path.exists(rootdir_codes)==False):
    os.chdir(rootdir_codes)


dir_in=[]
if(os.path.exists(rootdir_codes +'/opm-core')==False):
    os.chdir(rootdir_codes)
    a=os.system('git clone https://github.com/OPM/opm-core.git')
    if(a==0):
        dir_in=dir_in +['opm-core']    
else:
    os.chdir(rootdir_codes + '/opm-core')
    a=os.system('git pull')
    if(a==0):
         dir_in=dir_in +['opm-core']  

# define directories dependent of the installed
dir_prod=[]
if(os.path.exists(rootdir_codes + '/opm-polymer')==False):
    os.chdir(rootdir_codes)
    a=os.system('hg clone https://public.ict.sintef.no/opm/hg/opm-polymer')
    if(a==0):
        dir_prod=dir_prod +['opm-polymer']    
else:
    os.chdir(rootdir_codes + '/opm-polymer')
    os.system('hg pull -u')
    if(a==0):
        dir_prod=dir_prod +['opm-polymer']  

# define all directories
mydirs = dir_in + dir_prod

# use autoconf to initialize all projects
for mydir in mydirs:
   os.chdir(rootdir_codes + '/' + mydir)
   a=os.path.exists('configure')
   if(a==False):
       os.system('autoreconf -i')


# first configure the core dir which do not depend on the others
configure_dirs(rootdir_build,rootdir_codes,comopt,dir_in,False)   
            

# comple and make core dir
make_dirs(rootdir_build,comopt,dir_in,True)



# build all other projects which depend on opm-core
configure_dirs(rootdir_build,rootdir_codes,comopt,dir_prod,True)   
            
            
#compile                
# build all aother prodjects
make_dirs(rootdir_build,comopt,dir_prod,False)

