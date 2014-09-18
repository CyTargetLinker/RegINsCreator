package cytargetlinker.conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bridgedb.BridgeDb;
import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.graph.XGMMLWriter;
import cytargetlinker.conversion.utils.CommonAttributes;

/**
 * Downloads TransmiR file version 1.2
 * and converts it into a RegIN file to be used with 
 * CyTargetLinker
 * @author mkutmon
 *
 */
public class TransmiR {

	private static String transmiRUrl = "http://www.cuilab.cn/files/images/transmir/transmir_v1.2.txt";
	private static String bridgedbMappingFile = "/home/martina/Data/BridgeDb/Hs_Derby_20130701.bridge";
	private static String outputFile = "/home/martina/TransmiR-hsa-1.2.xgmml";
	
	public static void main(String[] args) {
		URL url;
		try {
			System.out.println("[INFO:]\tDownload TransmiR 1.2 file");
			url = new URL(transmiRUrl);
			URLConnection con = url.openConnection();
			BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()));
			reader.readLine();
			String line;
			List<String> list = new ArrayList<String>();
			while ((line = reader.readLine()) != null) {
				list.add(line);
			}
			reader.close();
			
			System.out.println("[INFO:]\tSet up identifier mapping.");
			File bridgeFile = new File(bridgedbMappingFile);
			Class.forName("org.bridgedb.rdb.IDMapperRdb");  
			IDMapper geneMapper = BridgeDb.connect("idmapper-pgdb:" + bridgeFile.getAbsolutePath());
			System.out.println("[INFO:]\tIdentifier mapping.database loaded");
			
			Graph graph = convert(list, geneMapper);
			writeGraph(graph, new File(outputFile));
			System.out.println("[INFO:]\tConversion is finished with " + geneNodes.size() + " TFs, " + miRNANodes.size() + " miRNAs and " + graph.getEdges().size() + " edges.");
			
		} catch (MalformedURLException e) {
			System.err.println("Can not read TransmiR file.\n" + e.getMessage());
		} catch (IOException e) {
			System.err.println("Can not read TransmiR file.\n" + e.getMessage());
		} catch (ClassNotFoundException e) {
			System.err.println("Can not load BridgeDb.\n" + e.getMessage());
		} catch (IDMapperException e) {
			System.err.println("Can not load BridgeDb databases.\n" + e.getMessage());
		}
	}
	
	private static Set<String> edges = new HashSet<String>();
	private static Map<String, Node> geneNodes = new HashMap<String, Node>();
	private static Map<String, Node> miRNANodes = new HashMap<String, Node>();
	
	private static Graph convert(List<String> interactions, IDMapper geneMapper) throws IDMapperException {
		Graph graph = new Graph();
		graph.setTitle("TransmiR v1.2");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "TransmiR v1.2");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "TF-miRNA interaction");
		
		for(String line : interactions) {
			String [] buffer = line.split("\t");
			String geneName = buffer[0];
			String entrez = buffer[1];
			String miRNA = buffer[3];
			if(buffer.length > 9) {
				String pubmed = buffer[8];
				String organism = buffer[9];
				if(organism.equals("human")) {
					Node source;
					if(geneNodes.containsKey(entrez)) {
						source = geneNodes.get(entrez);
					} else {
						source = addSourceNode(graph, entrez, geneMapper, geneName, organism);
						geneNodes.put(entrez, source);
					}
					
					
					Node target;
					if(miRNANodes.containsKey("hsa-" + miRNA)) {
						target = miRNANodes.get("hsa-" + miRNA);
					} else {
						target = graph.addNode("hsa-" + miRNA);
						target.appendAttribute("identifiers", "[" + "hsa-" + miRNA + "]");
						target.appendAttribute("label", "hsa-" + miRNA);
						target.appendAttribute("name", "hsa-" + miRNA);
						target.appendAttribute("biologicalType", "microRNA");
						miRNANodes.put("hsa-" + miRNA, target);
					}
					
					addEdge(source, target, graph, pubmed);
				}
			} else {
				System.out.println("No organism specified for interaction " + entrez + " -> " + miRNA);
			}
			
			
		}
		//gene	entrezid	tumor	mir	tumor_mir	mir_func	mir_disease	active	pmid	organism
		
		
		return graph;
	}
	
	private static Node addSourceNode(Graph graph, String entrez, IDMapper mapper, String geneName, String organism) throws IDMapperException {
		Node tf = graph.addNode(entrez);
		String identifiers = "[" + entrez;
		Xref x = new Xref(entrez, DataSource.getBySystemCode("L"));
		Set<Xref> res = mapper.mapID(x, DataSource.getBySystemCode("En"));
		for(Xref xref : res) {
			identifiers = identifiers + "," + xref.getId();
		}
		Set<Xref> res2 = mapper.mapID(x, DataSource.getBySystemCode("S"));
		for(Xref xref : res2) {
			identifiers = identifiers + "," + xref.getId();
		}
		identifiers = identifiers + "]";
		tf.appendAttribute("identifiers", identifiers);
		tf.appendAttribute("geneName", geneName);
		tf.appendAttribute("label", geneName);
		tf.appendAttribute("name", geneName);
		tf.appendAttribute("organism", organism);
		tf.appendAttribute("entrez", entrez);
		tf.appendAttribute("biologicalType", "transcriptionFactor");
		
		return tf;
	}
	
	private static void addEdge(Node source, Node target, Graph graph, String pubmed) {
		String id = source.getId() + "-" + target.getId();
		if(!edges.contains(id)) {
			Edge e = graph.addEdge(id, source, target);
			e.setAttribute("datasource", "TransmiR v1.2");
			e.setAttribute("interactionType", "TF-miRNA interaction");
			e.setAttribute("pubmed", pubmed);
		}
	}

	private static void writeGraph(Graph graph, File output) throws IOException {
		PrintWriter po = new PrintWriter(output);
		XGMMLWriter.write(graph, po);
		po.close();
	}
}
