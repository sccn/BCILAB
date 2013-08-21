/**
 *  To send a message to the path /target at destination d,
 *  int osc_send(d, struct('path', '/target', 'typetags', 'fff', 'data', {{1,1,2}}));
 *
 *  To send a bundle of messages, use a cell array of structures:
 *  There is currently no support for timestamps because I don't want to deal with leap seconds (ugh!).
 *  Nested bundles are not permitted.  Maximum of 256 nested messages in a bundle.
 *  int osc_send(d, { struct('path', '/target1', ...), struct(...), ... }) 
 *
 *  osc_send returns logical scalar 1 if it was successful (but this does not assure delivery due to UDP)
 *  otherwise returns 0 and should print an error (e.g. "Connection Denied")
**/

#include "mex.h"
#include "lo/lo.h"
#include "math.h"

#define EPOCH_DATENUM 693962.0

// Global buffers are faster but not thread-safe

char strbuf[4096];
char tphbuf[4096];
char dbgstr[256];
char msg_adr[256];
char msg_adri[256][256];
lo_message bdl_msgi[256];
char msg_tph[256];

void make_typehint(char* tph, int m, const mxArray *data) {

  int n, i;
  mxClassID t;

  if(data == NULL) { // allow for no data case
    tph[0] = '\0';
    return;   
  }
  
  if(! mxIsCell(data)) {
    mxErrMsgTxt("Expecting a cell array for data!");
  }
  n = mxGetNumberOfElements(data);

  for(i = 0; i < n && i < (m - 1); i++) {
    t = mxGetClassID(mxGetCell(data, i));

    switch(t) {
    case mxCHAR_CLASS:
      tph[i] = 's';
      break;
      
    case mxLOGICAL_CLASS:
      tph[i] = 'L';
      break;

    case mxDOUBLE_CLASS:
      tph[i] = 'd';
      break;

    case mxSINGLE_CLASS:
      tph[i] = 'f';
      break;

    case mxINT8_CLASS:
      tph[i] = 'c';
      break;

    case mxUINT8_CLASS:
      tph[i] = 'c';
      break;

    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:
      tph[i] = 'i';
      break;

    case mxINT64_CLASS:
    case mxUINT64_CLASS:
      tph[i] = 'h';
      break;

    default:
      mexErrMsgTxt("Sorry, the data contains a type that cannot be handled.");
      return;
    }
    
  }

  // terminate the type string...
  tph[i] = '\0';

}

void build_osc_message(lo_message msg, const char* typehint, const mxArray *data) {

  int n, m;
  int err;
  int i, j;
  int ti;

  const char* tph;
  
  lo_blob b;

  int d;

  if(data != NULL && ! mxIsCell(data)) {
    mexErrMsgTxt("Message data must be a cell array.");
    return;
  }

  if(data == NULL) {
    m = 0;
  } else {
    m = mxGetNumberOfElements(data);
  }

  if(! typehint) {
    make_typehint(tphbuf, 4096, data);
    tph = tphbuf;
  } else {
    tph = typehint;
  }

  n = strlen(tph);

  i = 0;
  ti = 0;

  d = 0;
  
  while(ti < n) {
    switch(tph[ti]) {
    case 'i':
      lo_message_add_int32(msg, (int)mxGetScalar(mxGetCell(data, i)));
      i++; ti++;
      break;
      
    case 'h':
      lo_message_add_int64(msg, (int64_t)mxGetScalar(mxGetCell(data, i)));
      i++; ti++;
      break;
      
    case 'f':
      lo_message_add_float(msg, (float)mxGetScalar(mxGetCell(data, i)));
      i++; ti++;
      break;
      
    case 'd':
      if(mxIsInf(mxGetScalar(mxGetCell(data, i)))) {
	lo_message_add_infinitum(msg);
      } else if(mxIsNaN(mxGetScalar(mxGetCell(data, i)))) {
	lo_message_add_nil(msg);
      } else {
	lo_message_add_double(msg, (double)mxGetScalar(mxGetCell(data, i)));
      }
      i++; ti++; break;
      
    case 's':
    case 'S':
    case 'b':
    case 'c':
      err = mxGetString(mxGetCell(data, i), strbuf, 4096);
      if(err != 0) {
	mexErrMsgTxt("Error reading string data in message contents (wrong typehint?)");
	return;
      }
      switch(tph[ti]) {
      case 's':
	lo_message_add_string(msg, strbuf);
	i++; ti++; break;
      case 'S':
	lo_message_add_symbol(msg, strbuf);
	i++; ti++; break;
      case 'b':
	mexErrMsgTxt("No blob support yet!"); return;
	lo_message_add_blob(msg, strbuf);
	i++; ti++; break;
      case 'c':
	lo_message_add_char(msg, strbuf[0]);
	i++; ti++; break;
      }
      break;
      
    case 'L': // "Logical" type (non-standard)
      if(mxIsLogicalScalar(mxGetCell(data, i)) && mxIsLogicalScalarTrue(mxGetCell(data, i)) || mxGetScalar(mxGetCell(data, i)) != 0) {
	lo_message_add_true(msg);
      } else {
	lo_message_add_false(msg);
      }
      i++; ti++;
      break;
      
    case 'T': // this group are the 'implicit' types, i.e. have zero data length
      lo_message_add_true(msg);
      ti++;
      break;
      
    case 'F':
      lo_message_add_false(msg);
      ti++;
      break;
      
    case 'N':
      lo_message_add_nil(msg);
      ti++;
      break;
      
    case 'I':
      lo_message_add_infinitum(msg);
      ti++;
      break;
      
    case '[':
    case ']':
      mexErrMsgTxt("No support for arrays!");
      return;

    case 'r':
      mexErrMsgTxt("No support for RGBA color type yet!");
      return;
    case 'm':
      mexErrMsgTxt("No support for MIDI type yet!");
      return;
      
    default:
      sprintf("typehint: %s at %d\n", tph, ti);
      mexErrMsgTxt("Unknown typehint.");
      return;
    }
    
  }
  
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  lo_address d;

  int i, j, n, m;
  int err;
  int mlen;

  const mxArray* p;
  const mxArray* q;
  const mxArray* r;
  
  lo_message msg;  // message pointer
  lo_bundle  bdl;  // bundle pointer if needed

  // no timestamp stuff for now...
  //lo_timetag ts;   // timestamp structure
  //double tsd;      // timestamp-as-double (in 'datenum' MATLAB format)

  if (nrhs != 2) {
    mexErrMsgTxt("Expecting two input arguments");
    return;
  }
  //if (nlhs != 1) {
  //  mexErrMsgTxt("Expecting one output arguments.");
  //  return;
  //}

  if(! mxIsChar(prhs[0])) {
    mexErrMsgTxt("Expecting a character array in the first argument.");
    return;
  }

  if(! (mxIsCell(prhs[1]) || mxIsStruct(prhs[1]))) {
    mexErrMsgTxt("Expecting a structure or a cell array in the second argument.");
  }

  mlen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  if(mlen > 32) {
    mexErrMsgTxt("Input string is too long.");
    return;
  }
  err = mxGetString(prhs[0], strbuf, mlen);

  if(err != 0) {
    mexErrMsgTxt("Error reading osc_address string.");
    return;
  }

  err = sscanf(strbuf, "osc_address:%x", &d);

  if(err < 0  || (! d)) {
    mexErrMsgTxt("Error scanning osc_address string.");
    return;
  }


  p = prhs[1];

  /*
  n = mxGetNumberOfElements(p);
  if(n == 0) {
    mexErrMsgTxt("Message structure must contain at least an address.");
  }

  err = mxGetString(mxGetCell(p, 0), msg_adr, 255);
  */
  
  if(mxIsCell(p)) { // its a bundle, i.e. list of structures

    /*  sorry no timestamp support for now...
    if(n == 3) {
      mexErrMsgTxt("Bundle timestamps are currently not supported.");
      tsd = mxGetScalar(mxGetCell(prhs[1], 1));
      ts.sec = (uint32_t)floor(tsd - EPOCH_DATENUM); // MATLAB timestamps are relative to Jan 1 0000, shift them to NTP standard.
      ts.frac = (uint32_t)((double)fmod(tsd, 1.0) * (4294967296.0)); // hmmm...
      bdl = lo_bundle_new(ts);

      sprintf(dbgstr, "bundle timestamp: %d, %f", (uint32_t)floor(tsd - EPOCH_DATENUM), fmod(tsd, 1.0));
      mexWarnMsgTxt(dbgstr);
    } else {
    }
    */
    
    bdl = lo_bundle_new(LO_TT_IMMEDIATE);

    n = mxGetNumberOfElements(p);

    for(i = 0; i < n && i < 256; i++) {
      r = mxGetCell(p, i);
      if(! mxIsStruct(r)) {
	mexErrMsgTxt("Bundled messages must be structures.");
	return;
      }

      bdl_msgi[i] = lo_message_new();
      msg_adri[i][0] = '\0';
      err = mxGetString(mxGetField(r, 0, "path"), msg_adri[i], 255);

      if(err != 0 || msg_adri[0] == NULL) {
	mexErrMsgTxt("Error reading address of bundled message.");
	for(i = 0; i < n; i++) {
	  lo_message_free(bdl_msgi[i]);
	}
	lo_bundle_free(bdl);
	return;
      }

      q = mxGetField(r, 0, "tt");
      msg_tph[0] = '\0';
      if(q != NULL) {
	  err = mxGetString(q, msg_tph, 255);
	  if(err != 0) {
	    mexErrMsgTxt("Error reading type hint of message.");
	  }
      }

      q = mxGetField(r, 0, "data");

      if(q || strlen(msg_tph)) {
	build_osc_message(bdl_msgi[i], strlen(msg_tph) ? msg_tph : NULL, q);
      }
      
      lo_bundle_add_message(bdl, msg_adri[i], bdl_msgi[i]);
    }
    // send bundle
    err = lo_send_bundle(d, bdl);
    if(err < 1) {
      mexWarnMsgTxt(lo_address_errstr(d));
      plhs[0] = mxCreateLogicalScalar(0);
    } else {
      plhs[0] = mxCreateLogicalScalar(1);
    }

    for(i = 0; i < n; i++) {
      lo_message_free(bdl_msgi[i]);
      msg_adri[i][0] = '\0'; // make sure these don't get reused
    }
    lo_bundle_free(bdl);

    return;
  } else {
    // its a normal message... (a struct)
    // parse cell and create message...

    msg = lo_message_new();

    msg_adr[0] = '\0';
    err = mxGetString(mxGetField(p, 0, "path"), msg_adr, 255);

    if(err != 0 || msg_adr == NULL) {
      mexErrMsgTxt("Error reading address of message.");
      lo_message_free(msg);
      return;
    }

    q = mxGetField(p, 0, "tt");
    msg_tph[0] = '\0';
    if(q != NULL) {
      err = mxGetString(q, msg_tph, 255);
      if(err != 0) {
	mexErrMsgTxt("Error reading type hint of message.");
      }
    }
    
    q = mxGetField(p, 0, "data");

    if(msg_adr && (q || strlen(msg_tph))) {
      build_osc_message(msg, strlen(msg_tph) ? msg_tph : NULL, q);
      err = lo_send_message(d, msg_adr, msg);

      if(err < 1) {
	mexWarnMsgTxt(lo_address_errstr(d));
	plhs[0] = mxCreateLogicalScalar(0);
	return;
      } else {
	plhs[0] = mxCreateLogicalScalar(1);
	return;
      }
    } else {
      mexErrMsgTxt("Error, no message to send!");
      return;
    }

    lo_message_free(msg);
  }
  
}

