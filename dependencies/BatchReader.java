import java.io.*;

class BatchReader {
    private InputStream in;
    public BatchReader(InputStream in) { this.in = in; }

    public byte[] readBatch(int num) throws Exception {
        byte[] buf = new byte[num];
        in.read(buf,0,num);
        return buf;
    }
}
