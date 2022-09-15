
'''
Read xcm files with mdefine commands and convert to slang functions.

Call as 'python parse_mdefine.py input_filename output_filename'
The output filename is optional, defaulting to replacing the input file suffix with '.sl'.

The list of additive fit functions may need extending for models outside than the base XSPEC set.
Functions defined and then used later in the same file are added automatically.

Some mathematical functions which may be used in mdefine statements are not included as ISIS standard functions; they must be defined before the output of your mdefine model is loaded.
The 'mdefine_specials' module includes all of these functions apart from 'dim'.

e.g.

ISIS> require("mdefine_specials");
ISIS> require("converted_xcm_file.sl");
ISIS> fit_fun("my_mdefine_function");
etc.

'''

import re
import numpy as np

## List of mathematical (not fitting) functions available to mdefine:
special_functions  = ['exp', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'sqrt', 'abs', \
    'asin', 'acos', 'atan' , 'asinh', 'acosh', 'atanh', 'sind', 'cosd', 'tand', \
        'heaviside', 'boxcar', 'sign', 'mean', 'atan2', 'erf', 'erfc', 'log', 'log10', \
            'gamma', 'min',' max', 'smin', 'smax']
## List of intrinsic XSPEC additive fit functions:
additive_functions = ['agauss', 'c6vmekl', 'eqpair', 'nei', 'rnei', 'vraymond', 'agnsed', 'carbatm', 'eqtherm', 'nlapec', 'sedov', 'vrnei', 'agnslim', 'cemekl', 'equil', 'npshock', 'sirf', 'vsedov', 'apec', 'cevmkl', 'expdec', 'nsa', 'slimbh', 'vtapec', 'bapec', 'cflow', 'ezdiskbb', 'nsagrav', 'smaug', 'vvapec', 'bbody', 'compLS', 'gadem', 'nsatmos', 'snapec', 'vvgnei', 'bbodyrad', 'compPS', 'gaussian', 'nsmax', 'srcut', 'vvnei', 'bexrav', 'compST', 'gnei', 'nsmaxg', 'sresc', 'vvnpshock', 'bexriv', 'compTT', 'grad', 'nsx', 'ssa', 'vvpshock', 'bkn2pow', 'compbb', 'grbcomp', 'nteea', 'step', 'vvrnei', 'bknpower', 'compmag', 'grbjet', 'nthComp', 'tapec', 'vvsedov', 'bmc', 'comptb', 'grbm', 'optxagn', 'vapec', 'vvtapec', 'bremss', 'compth', 'hatm', 'optxagnf', 'vbremss', 'vvwdem', 'brnei', 'cph', 'jet', 'pegpwrlw', 'vcph', 'vwdem', 'btapec', 'cplinear', 'kerrbb', 'pexmon', 'vequil', 'wdem', 'bvapec', 'cutoffpl', 'kerrd', 'pexrav', 'vgadem', 'zagauss', 'bvrnei', 'disk', 'kerrdisk', 'pexriv', 'vgnei', 'zbbody', 'bvtapec', 'diskbb', 'kyrline', 'plcabs', 'vmcflow', 'zbknpower', 'bvvapec', 'diskir', 'laor', 'posm', 'vmeka', 'zbremss', 'bvvrnei', 'diskline', 'laor2', 'powerlaw', 'vmekal', 'zcutoffpl', 'bvvtapec', 'diskm', 'logpar', 'pshock', 'vnei', 'zgauss', 'bwcycl', 'disko', 'lorentz', 'qsosed', 'vnpshock', 'zkerrbb', 'c6mekl', 'diskpbb', 'meka', 'raymond', 'voigt', 'zlogpar', 'c6pmekl', 'diskpn', 'mekal', 'redge', 'vpshock', 'zpowerlw', 'c6pvmkl', 'eplogpar', 'mkcflow', 'refsch']

def interpret_line(line,out):
    '''
    Convert an XSPEC mdefine model to ISIS/S-Lang syntax
    '''
    
    ### Separate line into components:

    ## Function name
    func_name = line.split(' ')[1]
    ## Function expression
    func_expr = ' '.join(line.split(':')[0].split(' ')[2:])
    ## Model type (add, mul, con); default to add if not given
    if len(line.split(':'))==1:
        mtype = 'add'
    else:
        mtype = line.split(':')[-1].strip()
    
    if mtype == 'add':
        additive_functions.append(func_name)
    
    ### Process XSPEC function expression:
    
    ## Replace '**' for exponentiation with '^':
    re.sub(r'\*\s*\*','^',func_expr)
    
    ## Correct logarithm bases:
    re.sub(r'\blog\b\s*\(','log10(',func_expr)
    re.sub(r'\bln\b\s*\(' ,'log('  ,func_expr)
    
    ## Change vector min/max syntax
    
    ## For binary min/max, put the two arguments in an array:
    for sub_func_name in ['min','max']:
            match = True
            new_func_expr = ''
        
            ##Find instance of the function:
            match = re.search(r'\b'+sub_func_name+r'\b', func_expr)
            while match:         
                
                # Move this segment to new string
                new_func_expr += func_expr[:match.start()]+sub_func_name+'(['
                func_expr = func_expr[match.end():]
                # Find string of the function arguments:                                
                func_expr = func_expr[re.search(r'\(',func_expr).end():]
                                                                                   
                i=0
                bracket_level = 1
                while bracket_level > 0:
                    if func_expr[i]==r'(':
                        bracket_level+=1
                    elif func_expr[i]==r')':
                        bracket_level-=1
                    i+=1
                func_expr=func_expr[:i-1]+'])'+func_expr[i:]
                                
                ##Find next instance of the function:
                match = re.search(r'\b'+sub_func_name+r'\b', func_expr)
                                
            func_expr = new_func_expr + func_expr
    
    re.sub(r'\bsmin\b\s*\(','min(',func_expr)
    re.sub(r'\bsmax\b\s*\(','max(',func_expr)
    
    ## Extract names of existing functions:
    before_bracket=[re.split('\W', _)[-1] for _ in func_expr.split('(') ]
    sub_func_names = [ _ for _ in before_bracket if _!='']
    
    ## Extract parameters (anything else that is not the energy specifier):
    words = [_ for _ in re.split('\W',func_expr) if _!='']
    words = [_ for _ in words if not re.match('[0-9]',_[0]) ]
    pars,inds = np.unique([ _ for _ in words if _ not in sub_func_names and _!='e' and _!='E' ], return_index=True)
    # Put parameters in order used in original function:
    pars     = list(pars[np.argsort(inds)])
    
    ### Make S-Lang expression by modifying XSPEC syntax:
        
    ## Address subfunctions properly and give wavelengths to subfunctions:
    for sub_func_name in np.unique(sub_func_names):
        if sub_func_name not in special_functions:
            match = True
            new_func_expr = ''
        
            ##Find instance of the function:
            match = re.search(r'\b'+sub_func_name+r'\b', func_expr)
            while match:         
                
                # Move this segment to new string
                new_func_expr += func_expr[:match.start()]+'eval_fun2(&'+sub_func_name+',lo,hi, ['
                func_expr = func_expr[match.end():]
                # Find string of the function arguments:
                
                # Add the normalisation for additive models:
                if sub_func_name in additive_functions:
                    new_func_expr += '1, '
                
                func_expr = func_expr[re.search(r'\(',func_expr).end():]
                                                                                   
                i=0
                bracket_level = 1
                while bracket_level > 0:
                    if func_expr[i]==r'(':
                        bracket_level+=1
                    elif func_expr[i]==r')':
                        bracket_level-=1
                    i+=1
                func_expr=func_expr[:i-1] +'])'+ func_expr[i:]
                
                ##Find next instance of the function:
                match = re.search(r'\b'+sub_func_name+r'\b', func_expr)
                
                
            func_expr = new_func_expr + func_expr
    
    '''
    Dealing with finding bin-integral quantities properly might need more care in the following two sections.
    '''
        
    ## Convert energy to ISIS expression (central energy):
    func_expr = re.sub(r'\b(e|E)\b', '(6.19920995*(lo+hi)/lo/hi)', func_expr)
    #(central wavelength):
    #func_expr = re.sub(r'\b(e|E)\b', '(24.7968398/(lo+hi))', func_expr)
    
    ## Allow a normalisation for additive models:
    if mtype == 'add':
        pars = ['norm']+pars
        func_expr = '( '+func_expr+' )*norm' # May need adjusting for units
    
    
    ### Make a compliant .sl file:
    
    ## Divide function expression if lines are too long:
    if len(func_expr)>72:
        func_expr = ' +\n        '.join(func_expr.split('+'))
        func_expr = ' *\n        '.join(func_expr.split('*'))
    
    ## Declare all parameters as variable names:
    declare_variables = 'variable '+', '.join(pars)
    if len(declare_variables)>72:
        declare_variables = 'variable '+', \n        '.join(pars)
    ## List of parameters for output:
    par_list = '["'+'","'.join(pars)+'"]'
    # Use separate lines if too long:
    if len(par_list)>72:
        par_list =  '["'+'",\n        "'.join(pars)+'"]'
    
    ### Write to file:
    out.write('define '+func_name+'_fit(lo,hi,par)\n{\n')
    out.write('    '+declare_variables+';\n')
    #out.write('    print("'+func_name+'");\n')
    for par,i in zip(pars,range(len(pars))):
        out.write('    '+par+' = par['+str(i)+'];\n')
    out.write('\n    return '+func_expr+';\n};\n\n')
    out.write('add_slang_function("'+func_name+'", '+par_list+');\n\n')
    
def convert_mdefine_file(input_xcm_file, output_sl_file):
    
    out = open(output_sl_file,'w')
    out.write('\n%%% Automatically translated by parse_mdefine (D. J. K. Buisson) %%%\n\n')
    with open(input_xcm_file,'r') as f:
        for line in f:
            if len(line.split())<1:
                pass
            elif line.split()[0] == '#':
                out.write('%'+line)
            elif line.split()[0] == 'mdefine':
                interpret_line(line,out)
            else:
                print('Failed to parse line: "'+line+'"')
    out.close()
    
if __name__ == "__main__":
    from os import sys
    ## Get filenames
    input_xcm_file = sys.argv[1]
    try:
        output_sl_file = sys.argv[2]
    except:
        ## Make output name if not given
        output_sl_file = '.'.join(input_xcm_file.split('.')[:-1])+'.sl'
    
    ## Make output
    convert_mdefine_file(input_xcm_file, output_sl_file)
