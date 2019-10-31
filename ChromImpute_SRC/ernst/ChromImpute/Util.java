package ernst.ChromImpute;

import java.io.*;
import java.util.*;
import java.util.zip.*;

public class Util
{
    /**
     * Returns a buffered reader. If szFile ends in a ".gz" tries to open it as a gzip file
     * otherwise tries to open it as a normal file.
     */
    public static BufferedReader getBufferedReader(String szFile) throws IOException
    {
        BufferedReader br;
	if (szFile.endsWith(".gz"))
	{
		br =new BufferedReader(new InputStreamReader(
							     new GZIPInputStream(new FileInputStream(szFile))));
        }
	else
	{
		br = new BufferedReader(new FileReader(szFile));
        }
	return br;
    }

}
