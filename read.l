%{
#include "read.tab.h"

int line_number = 1;
%}
%%
!.*\n { line_number++; }
\r?\n { line_number++; }
$DATA ;
$END ;
[ \t]+ ;
[+-]?[0-9]+[.][0-9]+ { yylval.f = atof(yytext); return FLOAT; }
[+-]?[0-9]+[.][0-9]+E[+-]?[0-9]+ { yylval.f = atof(yytext); return FLOAT; }
[+-]?[0-9]+ { yylval.i = atoi(yytext); return INTEGER; }
[A-Z][A-Z]+ { yylval.s = strdup(yytext); return ELEMENT; }
[A-Z] { yylval.i = yytext[0]; return SHELL; }
%%
int yywrap()
{
	return 1;
}

void yyerror(const char *s)
{
    fprintf(stderr, "error in line %d\n", line_number);
    exit(EXIT_FAILURE);
}
