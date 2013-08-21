import java.io.*;

class ChunkReader {
    private DataInput in;
    public ChunkReader(DataInput in) { this.in = in; }

    public byte[] readFully(int len) throws Exception {
        byte[] buf = new byte[len];
        in.readFully(buf);
        return buf;
    }
}