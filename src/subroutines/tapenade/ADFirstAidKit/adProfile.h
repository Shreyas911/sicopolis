#ifndef ADPROFILE_LOADED
#define ADPROFILE_LOADED

void adProfileAdj_endReverse(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_beginReverse(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_SNPWrite(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_SNPRead(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_beginAdvance(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_endAdvance(char* callname, char* filename, unsigned int lineno);
void adProfileAdj_turn(char* callname, char* filename);
void adProfileAdj_startRepeat() ;
void adProfileAdj_resetRepeat() ;
void adProfileAdj_endRepeat() ;
void adProfileAdj_showProfiles() ;


#endif // ADPROFILE_LOADED
