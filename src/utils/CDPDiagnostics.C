#include "CDPDiagnostics.h"
thread_local CDPDiagnostics::Context * CDPDiagnostics::current=nullptr;
thread_local double CDPDiagnostics::excluded_microseconds=0;
