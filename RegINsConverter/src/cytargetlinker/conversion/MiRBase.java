package cytargetlinker.conversion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.graph.XGMMLWriter;
import cytargetlinker.conversion.utils.CommonAttributes;

public class MiRBase {

	private static String mirbaseURL = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3";
	private static String outputFile = "/home/martina/mirbase-hsa-v21.xgmml";
	private static String networkName = "miRBase v 21";
	// you can use a biomart ensembl file to add links from immature to ensembl
	// use this line to specify a biomart file
//	private static String biomartFile = "/home/martina/Downloads/ensembl-mirna.txt";
	// use this line if you don't want to add a mapping 
	private static String biomartFile = "";
	
	public static void main(String[] args) throws Exception {
		
		MiRBase mirbase = new MiRBase();
		URL url = new URL(mirbaseURL);
		Map<String, Set<String>> ids = new HashMap<String, Set<String>>();
		if(!biomartFile.equals("")) {
			ids = readBiomartFile();
		}
		List<String> list = mirbase.readFile(url);
		Graph graph = mirbase.convertGraph(list, ids);
		mirbase.printGraph(graph);
	}
	
	

	private static Map<String, Set<String>> readBiomartFile() throws IOException {
		Map<String, Set<String>> map = new HashMap<String, Set<String>>();
		File f = new File(biomartFile);
		BufferedReader reader = new BufferedReader(new FileReader(f));
		String line;
		while((line = reader.readLine()) != null) {
			String [] buffer = line.split(",");
			if(buffer.length > 1) {
				String mirna = buffer[1];
				String ensembl = buffer[0];
				if(map.containsKey(mirna)) {
					if(!map.get(mirna).contains(ensembl)) {
						map.get(mirna).add(ensembl);
					}
				} else {
					Set<String> set = new HashSet<String>();
					set.add(ensembl);
					map.put(mirna, set);
				}
			}
		}
		reader.close();
		return map;
	}



	private Map<String, Node> genes;
	private Map<String, Node> mirnas;
	private Map<String, Edge> edges;
	
	public MiRBase() {
		genes = new HashMap<String, Node>();
		mirnas = new HashMap<String, Node>();
		edges = new HashMap<String, Edge>();
	}
	
	private void printGraph(Graph graph) throws Exception {
		PrintWriter po = new PrintWriter(new File(outputFile));
		XGMMLWriter.write(graph, po);
		po.close();
	}
	
	public List<String> readFile(URL url) throws Exception {
		System.out.println("[INFO]\treading mirBase file");
		URLConnection con = url.openConnection();
		BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()));
		reader.readLine();
		String line;
		List<String> list = new ArrayList<String>();
		while ((line = reader.readLine()) != null) {
			if(!line.startsWith("#")) {
				list.add(line);
			}
		}
		reader.close();
		System.out.println("[INFO]\tmirBase file has been read");
		return list;
	}

	public Graph convertGraph(List<String> list, Map<String, Set<String>> ids) {
		Graph graph = new Graph();
		graph.setTitle(networkName);
		graph.setAttribute(CommonAttributes.DATABASE.getName(), networkName);
		graph.setAttribute(CommonAttributes.TYPE.getName(), "primary transcript - miRNA interaction");
		for(String str : list) {
			String [] buffer = str.split("\t");
			if(buffer[2].equals("miRNA_primary_transcript")) {
				addSourceNode(buffer[8], graph, ids);		
			} else if(buffer[2].equals("miRNA")) {
				addTargetNode(buffer[8], graph);
			} 
		}
		System.out.println(networkName + " is converted to RegIN\n" + outputFile);
		System.out.println("Number of genes (primary transcripts): " + genes.size());
		System.out.println("Number of miRNAs: " + mirnas.size());
		System.out.println("Number of interactions: " + edges.size());
		return graph;
	}
	
	private void addTargetNode(String details, Graph graph) {
		String [] buffer2 = details.split(";");
		
		String mimat = buffer2[0].substring(3);
		String alias = buffer2[1].substring(6);
		String name = buffer2[2].substring(5);
		
		String gene = buffer2[3].substring(13);
		
		Node target;
		if(!mirnas.containsKey(mimat)) {
			target = graph.addNode(mimat);
			target.appendAttribute("miRBase id", mimat);
			target.appendAttribute("alias", alias);
			target.appendAttribute("name", name);
			target.appendAttribute("label", name);
			target.appendAttribute("biologicalType", "microRNA");
			target.appendAttribute("identifiers", "[" + mimat + "," + name + "]");
			mirnas.put(mimat, target);
		} else {
			target = mirnas.get(mimat);
		}
		
		Node geneNode = genes.get(gene);
		if(geneNode != null) {
			if(!edges.containsKey(geneNode.getId() + " - " + target.getId())) {
				Edge e = graph.addEdge(geneNode.getId() + " - " + target.getId(), geneNode, target);
				e.setAttribute("datasource", networkName);
				e.setAttribute("interactionType", "primary transcript - miRNA interaction");
				edges.put(geneNode.getId() + " - " + target.getId(), e);
			}
		} else {
			System.out.println("ERRROR!!!!");
		}
	}
	
	private void addSourceNode(String details, Graph graph, Map<String, Set<String>> ids) {
		String [] buffer2 = details.split(";");
		
		String mi = buffer2[0].substring(3);
		String alias = buffer2[1].substring(6);
		String name = buffer2[2].substring(5);
		
		if(!genes.containsKey(mi)) {
			String identifiers = "[" + mi + "," + name;
			if(ids.containsKey(mi)) {
				for(String ensembl : ids.get(mi)) {
					identifiers = identifiers + "," + ensembl;
				}
			}
			identifiers = identifiers + "]";
			System.out.println(identifiers);
			Node source = graph.addNode(mi);
			source.appendAttribute("miRBase id", mi);
			source.appendAttribute("alias", alias);
			source.appendAttribute("name", name);
			source.appendAttribute("label", name);
			source.appendAttribute("biologicalType", "gene");
			source.appendAttribute("identifiers", identifiers);
			genes.put(mi, source);
		}
	}
	
}
