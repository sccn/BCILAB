#include "bsb_common.h"
#include <time.h>

/*
 * [hDLL, hConn, hOPEN_DRIVER_ASYNC, hUSB_WRITE, hREAD_MULTIPLE_SWEEPS, hREAD_POINTER, hCLOSE_DRIVER_ASYNC, pBuffer, nbchan, srate, numSync, numTri, numEEG, numExG, last_ptr] = bsb_open(dllpath);
 * Backend function for BioSemi interface; do not use directly (instead via bs_open)
 * 
 * In:
 *   dllpath : path to the dll to be loaded (usually in the same dir as this file)
 *
 * Out:
 *   hDLL    : handle to the driver DLL
 *   hConn   : connection handle from the driver 
 *   hOPEN_DRIVER_ASYNC    : handle to a driver function
 *   hUSB_WRITE            : handle to a driver function
 *   hREAD_MULTIPLE_SWEEPS : handle to a driver function
 *   hREAD_POINTER         : handle to a driver function
 *   hCLOSE_DRIVER_ASYNC   : handle to a driver function
 *   pBuffer : pointer to the circular buffer (int*)
 *   nbchan  : number of channels per sample (consists of Sync,Trigger,EEG,ExG and Extra channels)
 *   srate   : sampling rate
 *   numSync  : number of synch channels among nbchan
 *   numTri  : number of trigger channels among nbchan
 *   numEEG  : number of EEG channels among nbchan
 *   numExG  : number of ExG channels among nbchan
 *   last_ptr : last buffer offset (to be used as last_ptr in bsb_read)
 */


/* length of the probing buffer in bytes */
#define BUFFER_BYTES (8*1024*1024)
#define MAX_WAITING_TIME 3.0

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) 
{
    void *hDLL;             /* handle to the DLL */
    char dllpath[16384];    /* path from where we load the dll (first arg) */
    int dllpath_len;
    /* function pointers */
    OPEN_DRIVER_ASYNC_t OPEN_DRIVER_ASYNC;
    USB_WRITE_t USB_WRITE;
    READ_MULTIPLE_SWEEPS_t READ_MULTIPLE_SWEEPS;
    READ_POINTER_t READ_POINTER;
    CLOSE_DRIVER_ASYNC_t CLOSE_DRIVER_ASYNC;
    /* comm data */
    void *hConn;            /* handle to the driver connection */
    char msg_enable[64] = {0};   /* 64 zero bytes */
    unsigned char msg_handshake[64] = {255,0};   /* 0xFF followed by 63 zero bytes */
    uint32_T  *ring_buffer;  /* pointer to ring buffer (once initialized) */
    /* ring buffer probing */
    uint32_T start_ptr;     /* pointer to the first sample that was recorded */
    uint32_T cur_ptr;       /* pointer past the current sample */
    clock_t start_time;     /* start time... */
    /* status info */
    uint32_T status;        /* content of the status flag */
    int32_T is_mk2;         /* whether the amp is a MK2 amplifier */
    int32_T speed_mode;     /* amplifier speed mode */
    int32_T srate;          /* sampling rate */
    int32_T nbsync;         /* number of synchronization channels */
    int32_T nbtrig;         /* number of trigger channels */
    int32_T nbeeg;          /* number of EEG channels */
    int32_T nbexg;          /* number of ExG channels */
    int32_T nbextra;        /* number of extra channels */
    int32_T nbchan;         /* total number of channels */
    int32_T battery_low;    /* whether the battery is low */
    /* result preparation */
    mwSize scalar[] = {1,1};
    int32_T last_ptr = 0; /* last pointer position (by default: 0) */
   
    /* check input/output arguments */   
    if (nrhs != 1)
        mexErrMsgTxt("One input argument required."); 
    if (nlhs != 15)
        mexErrMsgTxt("15 output arguments required."); 
    if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1) || (mxGetN(prhs[0]) <= 4))
        mexErrMsgTxt("First input argument must be a string."); 
    
    /* get DLL path */
    dllpath_len = mxGetNumberOfElements(prhs[0]);
    mxGetNChars_700(prhs[0], dllpath, dllpath_len+1);
    
    /* load the library */
    mexPrintf("Loading BioSemi driver dll...\n");
    hDLL = (void*)LOAD_LIBRARY(dllpath);
    if (!hDLL)
        mexErrMsgTxt("Could not load BioSemi driver DLL.");
    /* obtain DLL functions */
    OPEN_DRIVER_ASYNC = (OPEN_DRIVER_ASYNC_t)LOAD_FUNCTION(hDLL,"OPEN_DRIVER_ASYNC");
    if (!OPEN_DRIVER_ASYNC) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Did not find function OPEN_DRIVER_ASYNC() in the BisoSemi driver DLL.");
    }
    USB_WRITE = (USB_WRITE_t)LOAD_FUNCTION(hDLL,"USB_WRITE");
    if (!USB_WRITE) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Did not find function USB_WRITE() in the BisoSemi driver DLL.");
    }
    READ_MULTIPLE_SWEEPS = (READ_MULTIPLE_SWEEPS_t)LOAD_FUNCTION(hDLL,"READ_MULTIPLE_SWEEPS");
    if (!READ_MULTIPLE_SWEEPS) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Did not find function READ_MULTIPLE_SWEEPS() in the BisoSemi driver DLL.");
    }
    READ_POINTER = (READ_POINTER_t)LOAD_FUNCTION(hDLL,"READ_POINTER");
    if (!READ_POINTER) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Did not find function READ_POINTER() in the BisoSemi driver DLL.");
    }
    CLOSE_DRIVER_ASYNC = (CLOSE_DRIVER_ASYNC_t)LOAD_FUNCTION(hDLL,"CLOSE_DRIVER_ASYNC");
    if (!CLOSE_DRIVER_ASYNC) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Did not find function CLOSE_DRIVER_ASYNC() in the BisoSemi driver DLL.");
    }
    
    /* connect to driver */
    mexPrintf("Connecting to driver...\n");
    hConn = OPEN_DRIVER_ASYNC();
    if (!hConn) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Could not open connection to BioSemi driver.");
    }
    
    /* initialize USB2 interface */
    mexPrintf("Initializing USB interface...\n");
    if (!USB_WRITE(hConn,msg_enable)) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Could not initialize BioSemi USB2 interface.");
    }
    
    /* initialize the initial (probing) ring buffer */
    mexPrintf("Initializing ring buffer...\n");
    ring_buffer = (int*) malloc(BUFFER_BYTES);
    if (!ring_buffer) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);        
        mexErrMsgTxt("Could not allocate ring buffer (out of memory).");
    } 
    memset(ring_buffer,0,BUFFER_BYTES);
    
    /* begin acquisition */
    READ_MULTIPLE_SWEEPS(hConn,(char*)ring_buffer,BUFFER_BYTES);
    
    /* enable handshake */
    mexPrintf("Enabling handshake...\n");
    if (!USB_WRITE(hConn,msg_handshake)) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("Could not enable handshake with BioSemi USB2 interface.");
    }
    
    /* obtain buffer head pointer */
    mexPrintf("Querying buffer pointer...\n");
    if (!READ_POINTER(hConn,&start_ptr)) {
        USB_WRITE(hConn, msg_enable);
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("Can not obtain ring buffer pointer from BioSemi driver.");
    }
    
    /* check head pointer validity */
    if (start_ptr > BUFFER_BYTES) {
        USB_WRITE(hConn, msg_enable);
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("Buffer pointer returned by BioSemi driver is not in the valid range.");
    }
    
    /* read the first sample */
    mexPrintf("Waiting for data...\n");
    start_time = clock();
    while (1) {
        if (!READ_POINTER(hConn, &cur_ptr)) {
            USB_WRITE(hConn, msg_enable);
            CLOSE_DRIVER_ASYNC(hConn);
            FREE_LIBRARY(hDLL);
            free(ring_buffer);
            mexErrMsgTxt("Ring buffer handshake with BioSemi driver failed unexpectedly.");
        }
        if ( ((double)(clock() - start_time))/CLOCKS_PER_SEC > MAX_WAITING_TIME) {
            USB_WRITE(hConn, msg_enable);
            CLOSE_DRIVER_ASYNC(hConn);
            FREE_LIBRARY(hDLL);
            free(ring_buffer);
            if (cur_ptr - start_ptr < 8)
                mexErrMsgTxt("BioSemi driver does not transmit data. Is the box turned on?");
            else
                mexErrMsgTxt("Did not get a sync signal from BioSemi driver. Is the battery charged?");
        }
        if ((cur_ptr - start_ptr >= 8) && (ring_buffer[start_ptr/4] == 0xFFFFFF00))
            break;
    }
    
    /* read the trigger channel data ... */
    mexPrintf("Checking status...\n");
    status = ring_buffer[start_ptr/4+1] >> 8;
    mexPrintf("  status: %u\n",status);

    /* determine channel configuration */
    is_mk2 = ((status&(1<<23)) != 0);    
    mexPrintf("  MK2: %u\n",is_mk2);
    
    /* check speed mode */    
    speed_mode = ((status&(1<<17))>>17) + ((status&(1<<18))>>17) + ((status&(1<<19))>>17) + ((status&(1<<21))>>18);
    mexPrintf("  speed mode: %u\n",speed_mode);

    if (!is_mk2 && (speed_mode == 8)) {
        USB_WRITE(hConn, msg_enable);
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("BioSemi amplifier speed mode wheel must be between positions 0 and 7 (currently 8); recommended is 4.");
    }
    if (speed_mode > 9) {
        USB_WRITE(hConn, msg_enable);
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        if (is_mk2)
            mexErrMsgTxt("BioSemi amplifier speed mode wheel must be between positions 0 and 8 (currently 9 = Reserved); recommended is 4.");
        else
            mexErrMsgTxt("BioSemi amplifier speed mode wheel must be between positions 0 and 7 (currently 9 = Reserved); recommended is 4.");
    }
    
    /* determine sampling rate */
    switch (speed_mode & 3) {
        case 0: srate = 2048; break;
        case 1: srate = 4096; break;
        case 2: srate = 8192; break;
        case 3: srate = 16384; break;
    }
    mexPrintf("  sampling rate: %u\n",srate);

    if (is_mk2) {
        switch (speed_mode) {
            case 0: 
            case 1: 
            case 2: 
            case 3:
                nbeeg = 0; nbexg = 4*152; nbtrig = 1; nbsync = 1; nbchan = 610; break;
            case 4: nbeeg = 256; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 282; break;
            case 5: nbeeg = 128; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 154; break;
            case 6: nbeeg = 64; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 90; break;
            case 7: nbeeg = 32; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 314; break;
            case 8: nbeeg = 280; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 314; break;
        }
    } else {
        switch (speed_mode) {
            case 0: nbeeg = 256; nbexg = 0; nbtrig = 1; nbsync = 1; nbchan = 258; break;
            case 1: nbeeg = 128; nbexg = 0; nbtrig = 1; nbsync = 1; nbchan = 130; break;
            case 2: nbeeg = 64; nbexg = 0; nbtrig = 1; nbsync = 1; nbchan = 66; break;
            case 3: nbeeg = 32; nbexg = 0; nbtrig = 1; nbsync = 1; nbchan = 34; break;
            case 4: nbeeg = 232; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 258; break;
            case 5: nbeeg = 104; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 130; break;
            case 6: nbeeg = 40; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 66; break;
            case 7: nbeeg = 8; nbexg = 8; nbtrig = 1; nbsync = 1; nbchan = 34; break;
        }        
    }
    nbextra = nbchan - (nbeeg + nbexg + nbtrig + nbsync);
    mexPrintf("  channels: %u (%uxSync, %uxTrigger, %uxEEG, %uxExG, %uxExtra)\n",nbchan,nbsync,nbtrig,nbeeg,nbexg,nbextra);

    /* check for additional problematic conditions */
    battery_low = status & (1<<22);
    if (battery_low)
        mexPrintf("  battery: low\n");
    else
        mexPrintf("  battery: good\n");
    if (battery_low)
        mexWarnMsgIdAndTxt("BioSemi:BatteryLow","The BioSemi battery is low; amplifier will shut down within 30-60 minutes.");
     
    /* compute a proper buffer size (needs to be a multiple of 512, a multiple of nbchan, as well as ~32MB in size. */
    mexPrintf("Reallocating optimized ring buffer...\n");

    /* shutdown current connection */
    mexPrintf("Sending the enable message again...\n");
    if (!USB_WRITE(hConn, msg_enable)) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Error while disabling the handshake.");
    }
    mexPrintf("Closing the driver...\n");
    if (!CLOSE_DRIVER_ASYNC(hConn)) {
        FREE_LIBRARY(hDLL);        
        mexErrMsgTxt("Error while disconnecting.");
    } 
            
    /* reset to a new ring buffer */
    mexPrintf("Freeing the ring buffer...\n");
    free(ring_buffer);
    mexPrintf("Allocating a new ring buffer...\n");
    ring_buffer = malloc(BUFFER_SAMPLES*4*nbchan);
    if (!ring_buffer) {
        FREE_LIBRARY(hDLL);        
        mexErrMsgTxt("Could not reallocate ring buffer (out of memory?).");
    } 
    mexPrintf("Zeroing the ring buffer...\n");
    memset(ring_buffer,0,BUFFER_SAMPLES*4*nbchan);
    
    /* reconnect to driver */
    mexPrintf("Opening the driver...\n");
    hConn = OPEN_DRIVER_ASYNC();
    if (!hConn) {
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Could not open connection to BioSemi driver.");
    }
    /* reinitialize USB2 interface */
    mexPrintf("Reinitializing the USB interface...\n");
    if (!USB_WRITE(hConn,msg_enable)) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        mexErrMsgTxt("Could not initialize BioSemi USB2 interface.");
    }    
    /* begin acquisition */
    mexPrintf("Starting data acquisition...\n");
    READ_MULTIPLE_SWEEPS(hConn,(char*)ring_buffer,BUFFER_SAMPLES*4*nbchan);    
    /* enable handshake */
    mexPrintf("Enabling handshake...\n");
    if (!USB_WRITE(hConn,msg_handshake)) {
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("Could not reenable handshake with BioSemi USB2 interface.");
    }
    /* make sure that we can read the buffer pointer */
    mexPrintf("Attempt to read buffer pointer...\n");
    if (!READ_POINTER(hConn,&start_ptr)) {
        USB_WRITE(hConn, msg_enable);
        CLOSE_DRIVER_ASYNC(hConn);
        FREE_LIBRARY(hDLL);
        free(ring_buffer);
        mexErrMsgTxt("Can not obtain ring buffer pointer from BioSemi driver.");
    }
    
    mexPrintf("Verifying channel format...\n");
    start_time = clock();
    while (1) {
        if (!READ_POINTER(hConn, &cur_ptr)) {
            USB_WRITE(hConn, msg_enable);
            CLOSE_DRIVER_ASYNC(hConn);
            FREE_LIBRARY(hDLL);
            free(ring_buffer);
            mexErrMsgTxt("Ring buffer handshake with BioSemi driver failed unexpectedly.");
        }
        if ( ((double)(clock() - start_time))/CLOCKS_PER_SEC > MAX_WAITING_TIME) {
            USB_WRITE(hConn, msg_enable);
            CLOSE_DRIVER_ASYNC(hConn);
            FREE_LIBRARY(hDLL);
            free(ring_buffer);
            if (cur_ptr - start_ptr < 8)
                mexErrMsgTxt("BioSemi driver does not transmit data. Is the box turned on?");
            else
                mexErrMsgTxt("Did not get a sync signal from BioSemi driver. Is the battery charged?");
        }
        if ((cur_ptr - start_ptr >= 4*nbchan) && (ring_buffer[start_ptr/4] == 0xFFFFFF00)) {
            if (ring_buffer[start_ptr/4+nbchan] != 0xFFFFFF00) {
                USB_WRITE(hConn, msg_enable);
                CLOSE_DRIVER_ASYNC(hConn);
                FREE_LIBRARY(hDLL);
                free(ring_buffer);
                mexErrMsgTxt("Sync signal did not show up at the expected position.");
            } else {
                mexPrintf("Channel format is correct...\n");
                break;
            }
        }
    }    
        
    mexPrintf("Setting up MATLAB return variables...\n");
    /* convert to MATLAB return values... */
    /*   hDLL    : handle to the driver DLL (int64) */
    plhs[0] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[0]) = (uintptr_t)hDLL;
    /*   hConn   : connection handle from the driver (int64) */
    plhs[1] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[1]) = (uintptr_t)hConn;    
    /*   hOPEN_DRIVER_ASYNC    : handle to a driver function */
    plhs[2] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[2]) = (uintptr_t)OPEN_DRIVER_ASYNC;
    /*   hUSB_WRITE            : handle to a driver function */
    plhs[3] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[3]) = (uintptr_t)USB_WRITE;
    /*   hREAD_MULTIPLE_SWEEPS : handle to a driver function */
    plhs[4] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[4]) = (uintptr_t)READ_MULTIPLE_SWEEPS;
    /*   hREAD_POINTER         : handle to a driver function */
    plhs[5] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[5]) = (uintptr_t)READ_POINTER;
    /*   hCLOSE_DRIVER_ASYNC   : handle to a driver function */
    plhs[6] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[6]) = (uintptr_t)CLOSE_DRIVER_ASYNC;
    /*   pBuffer : pointer to the circular buffer (int64,int*) */
    plhs[7] = mxCreateNumericArray(2,scalar,PTR_CLASS,mxREAL); *(uintptr_t*)mxGetData(plhs[7]) = (uintptr_t)ring_buffer;
    /*   nbchan  : number of channels per sample */
    plhs[8] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[8]) = (int32_T)nbchan;
    /*   srate   : sampling rate */
    plhs[9] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[9]) = (int32_T)srate;
    /*   numSync : number of synch channels among nbchan */
    plhs[10] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[10]) = (int32_T)nbsync;
    /*   numTri  : number of trigger channels among nbchan */
    plhs[11] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[11]) = (int32_T)nbtrig;
    /*   numEEG  : number of EEG channels among nbchan */
    plhs[12] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[12]) = (int32_T)nbeeg;
    /*   numExG  : number of ExG channels among nbchan */
    plhs[13] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[13]) = (int32_T)nbexg;    
    /*   last_ptr : last buffer offset */
    plhs[14] = mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL); *(int32_T*)mxGetData(plhs[14]) = (int32_T)last_ptr;
}
