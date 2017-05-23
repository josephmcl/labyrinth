"""
A simple pybison parser program implementing a calculator
"""
from bison import BisonParser


class Parser(BisonParser):
    """
    Implements the calculator parser. Grammar rules are defined in the method
    docstrings. Scanner rules are in the 'lexscript' attribute.
    """
    # ----------------------------------------------------------------
    # lexer tokens - these must match those in your lex script (below)
    # ----------------------------------------------------------------
    tokens = ['CLASS','DEF','EXTENDS','IF','ELIF','ELSE','WHILE','RETURN',
              'IDENT', 'INT_LIT', 'TIMES', 'DIVIDE', 'POW',
              'LPAREN', 'RPAREN',
              'NEWLINE', 'QUIT']

    # ------------------------------
    # precedences
    # ------------------------------
    precedences = (
        ('left', ('MINUS', 'PLUS')),
        ('left', ('TIMES', 'DIVIDE')),
        ('left', ('NEG', )),
        ('right', ('POW', )),
        )

    # ------------------------------------------------------------------
    # override default read method with a version that prompts for input
    # ------------------------------------------------------------------
    def read(self, nbytes):
        try:
            return raw_input("> ") + "\n"
        except EOFError:
            return ''

    # ---------------------------------------------------------------
    # These methods are the python handlers for the bison targets.
    # (which get called by the bison code each time the corresponding
    # parse target is unambiguously reached)
    #
    # WARNING - don't touch the method docstrings unless you know what
    # you are doing - they are in bison rule syntax, and are passed
    # verbatim to bison to build the parser engine library.
    # ---------------------------------------------------------------

    # Declare the start target here (by name)
    start = "input"

    def on_input(self, target, option, names, values):
        """
        input :
              | input line
        """
        return

    def on_line(self, target, option, names, values):
        """
        line : NEWLINE
             | exp NEWLINE
        """
        if option == 1:
            print values[0]

    def on_exp(self, target, option, names, values):
        """
        exp : NUMBER
            | exp PLUS exp
            | exp MINUS exp
            | exp TIMES exp
            | exp DIVIDE exp
            | MINUS exp %prec NEG
            | exp POW exp
            | LPAREN exp RPAREN
        """
        #print "on_exp: got %s %s %s %s" % (target, option, names, values)
        if option == 0:
            return float(values[0])
        elif option == 1:
            return values[0] + values[2]
        elif option == 2:
            return values[0] - values[2]
        elif option == 3:
            return values[0] * values[2]
        elif option == 4:
            return values[0] / values[2]
        elif option == 5:
            return - values[1]
        elif option == 6:
            return values[0] ** values[2]
        elif option == 7:
            return values[1]

    # -----------------------------------------
    # raw lex script, verbatim here
    # -----------------------------------------
    lexscript = r"""
    %{
    //int yylineno = 0;
    #include <stdio.h>
    #include <string.h>
    #include "Python.h"
    #define YYSTYPE void *
    #include "tokens.h"
    extern void *py_parser;
    extern void (*py_input)(PyObject *parser, char *buf, int *result,
                            int max_size);
    #define returntoken(tok) \
            yylval = PyString_FromString(strdup(yytext)); return (tok);
    #define YY_INPUT(buf,result,max_size) { \
        (*py_input)(py_parser, buf, &result, max_size); \
    }
    %}
    %%
    %option yylineno

    %x BLOCK_COMMENT
    %%

    \n                  {}
    [ \t\r\f]               /* do nothing */

    "/*"                        {BEGIN(BLOCK_COMMENT);}
    <BLOCK_COMMENT>"*/"             {BEGIN(INITIAL);}
    <BLOCK_COMMENT>\n               {}
    <BLOCK_COMMENT>.                {}


    "\"\"\""([^']|'[^']|''[^'])*"\"\"\""    {yylval.sval = strdup(yytext); return STRING_LIT;}

    [/]{2}+.*               {}
    \"(\\[0bftnr\\]|[^"\n\\])*\"        {yylval.sval = strdup(yytext); return STRING_LIT;}
    \"(\\.|[^"]|\n)*\"          {return 21;}


    class                   {yylval.sval = strdup(yytext); return CLASS;}
    def                 {yylval.sval = strdup(yytext); return DEF;}
    extends                 {yylval.sval = strdup(yytext); return EXTENDS;}
    if                  {yylval.sval = strdup(yytext); yylval.sval = strdup(yytext); return IF;}
    elif                    {yylval.sval = strdup(yytext); return ELIF;}
    else                    {yylval.sval = strdup(yytext); return ELSE;}
    while                   {yylval.sval = strdup(yytext); return WHILE;}
    return                  {yylval.sval = strdup(yytext); return RETURN;}

    String                  {yylval.sval = strdup(yytext); return IDENT;}   
    Integer                 {yylval.sval = strdup(yytext); return IDENT;}
    Obj                 {yylval.sval = strdup(yytext); return IDENT;}
    Boolean                 {yylval.sval = strdup(yytext); return IDENT;}
    true                    {yylval.sval = strdup(yytext); return IDENT;}
    false                   {yylval.sval = strdup(yytext); return IDENT;}
    and                 {yylval.sval = strdup(yytext); return AND;}
    or                  {yylval.sval = strdup(yytext); return OR;}
    not                 {yylval.sval = strdup(yytext); return NOT;}
    Nothing                 {yylval.sval = strdup(yytext); return IDENT;}
    none                    {yylval.sval = strdup(yytext); return IDENT;}

    [_a-zA-Z][_a-zA-Z0-9]*          {yylval.sval = strdup(yytext); return IDENT;}
    [1-9][0-9]*             {yylval.sval = strdup(yytext); return INT_LIT;}
    0                   {yylval.sval = strdup(yytext); return INT_LIT;}

    "+"                 {yylval.sval = strdup(yytext); return '+';}
    "-"                 {yylval.sval = strdup(yytext); return '-';}
    "*"                 {yylval.sval = strdup(yytext); return '*';}
    "/"                 {yylval.sval = strdup(yytext); return '/';}

    "="{2}                  {yylval.sval = strdup(yytext); return EQUALS;}
    "<="                    {yylval.sval = strdup(yytext); return ATMOST;}
    "<"                 {yylval.sval = strdup(yytext); return '<';}
    ">="                    {yylval.sval = strdup(yytext); return ATLEAST;}
    ">"                 {yylval.sval = strdup(yytext); return '>';}

    "{"                 {yylval.sval = strdup(yytext); return '{';}
    "}"                 {yylval.sval = strdup(yytext); return '}';}
    "("                 {yylval.sval = strdup(yytext); return '(';}
    ")"                 {yylval.sval = strdup(yytext); return ')';}
    ","                 {yylval.sval = strdup(yytext); return ',';}
    "."                 {yylval.sval = strdup(yytext); return '.';}
    ";"                 {yylval.sval = strdup(yytext); return ';';}
    ":"                 {yylval.sval = strdup(yytext); return ':';}
    "="                 {yylval.sval = strdup(yytext); return '=';}

    <<EOF>>                 {return 0;}

    .                   {return 20;}
    """

if __name__ == '__main__':
    p = Parser()
p.run()