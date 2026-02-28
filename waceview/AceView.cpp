/*  $Id: AceView.cpp 130374 2008-06-10 11:05:01Z grichenk $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Denis Vakatov, 
 *          Modifed by Jean Thierry-Mieg: wrapper for an external C cgi 
 *
 * File Description:
 *   CGI/FastCGI application, serving as a wrapper for external C code constructing an html page
 *
 *   USAGE:  sample.cgi?message=Some+Message
 *
 *   NOTE:
 *     1) needs the externally compiled aceviewlib module and libraries
 *
 */

#include <ncbi_pch.hpp>
#include <cgi/cgiapp.hpp>
#include <cgi/cgictx.hpp>

// aceviewlib(query) export a complete html page, given the webbrowser query
extern "C" char *aceviewlib (void) ;

using namespace ncbi;

/////////////////////////////////////////////////////////////////////////////
//  CCgiSampleApplication::
//

class CCgiSampleApplication : public CCgiApplication
{
public:
    virtual void Init(void);
    virtual int  ProcessRequest(CCgiContext& ctx);
};


void CCgiSampleApplication::Init()
{
    // Standard CGI framework initialization
    CCgiApplication::Init();
}


int CCgiSampleApplication::ProcessRequest(CCgiContext& ctx)
{
    // Given "CGI context", get access to its "HTTP response" sub-object
    CCgiResponse&      response = ctx.GetResponse();

    // Call the legacy C code to create the HTML page and flush it
    response.out() <<  aceviewlib() ;
    response.Flush();

    _TRACE("Hello world") ;  // no idea where this is going to
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
//  MAIN
//

int main(int argc, const char* argv[])
{
    GetDiagContext().SetOldPostFormat(false); // Switch to the new log format
    int result = CCgiSampleApplication().AppMain(argc, argv, 0, eDS_Default);
    _TRACE("back to normal diags");
    return result;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

