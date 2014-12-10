package cytargetlinker.conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.XGMMLWriter;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.utils.CommonAttributes;

public class TabFileConverter {

	private static String inputFile = "resources/Correlation.csv";
	private static String outputFile = "resources/Correlation.xgmml";
	private static String networkName = "Correlation";
	
	private static int [] sourceIdCol = {0,1};
	private static int [] targetIdCol = {8};
	private static int [] sourceAttrCol = {1,2,3,4,5,6,7};
	private static int [] targetAttrCol = {};
	private static int [] edgeAttrCol = {};
	
	public static void main(String[] args) throws Exception {
		
		Graph graph = new Graph();
		graph.setTitle(networkName);
		graph.setAttribute(CommonAttributes.DATABASE.getName(), networkName);		
		
		File f = new File(inputFile);
		BufferedReader reader = new BufferedReader(new FileReader(f));
		String line;
		String [] header = reader.readLine().split("\t");
		while((line = reader.readLine()) != null) {
			String [] buffer = line.split("\t");
			String sourceId = buffer[sourceIdCol[0]];
			Node source = graph.addNode(sourceId);
			String sourceIdentifier = "[" + sourceId;
			for(int i = 1; i < sourceIdCol.length; i++) {
				sourceIdentifier = sourceIdentifier + "," + buffer[sourceIdCol[i]];
			}
			sourceIdentifier = sourceIdentifier + "]";
			source.appendAttribute("identifiers", sourceIdentifier);
			for(int i = 0; i < sourceAttrCol.length; i++) {
				source.appendAttribute(header[sourceAttrCol[i]], buffer[sourceAttrCol[i]]);
			}

			String id = buffer[targetIdCol[0]];
			Node target = graph.addNode(id);
			String identifier = "[" + id;
			for(int i = 1; i < targetIdCol.length; i++) {
				identifier = identifier + "," + buffer[targetIdCol[i]];
			}
			identifier = identifier + "]";
			target.appendAttribute("identifiers", identifier);
			for(int i = 0; i < targetAttrCol.length; i++) {
				target.appendAttribute(header[targetAttrCol[i]], buffer[targetAttrCol[i]]);
			}
		
			Edge e = graph.addEdge(sourceId + "-" + id, source, target);
			for(int i = 0; i < edgeAttrCol.length; i++) {
				e.appendAttribute(header[edgeAttrCol[i]], buffer[edgeAttrCol[i]]);
			}
			
		}
		 
		reader.close();
		
		writeGraph(graph,new File(outputFile));
	}
	
	private static void writeGraph(Graph graph, File output) throws IOException {
		PrintWriter po = new PrintWriter(output);
		XGMMLWriter.write(graph, po);
		po.close();
	}

}
