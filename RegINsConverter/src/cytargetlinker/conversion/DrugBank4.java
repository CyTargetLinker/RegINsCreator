package cytargetlinker.conversion;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import msk.drugbank4.DrugBankParser;
import msk.drugbank4.DrugModel;
import msk.drugbank4.TargetModel;

import org.bridgedb.BridgeDb;
import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;
import org.jdom.JDOMException;

import cytargetlinker.conversion.graph.Graph;
import cytargetlinker.conversion.graph.Graph.Edge;
import cytargetlinker.conversion.graph.Graph.Node;
import cytargetlinker.conversion.graph.XGMMLWriter;
import cytargetlinker.conversion.utils.CommonAttributes;

public class DrugBank4 {

	public static void main(String[] args) {
		File drugBankXml = new File("/home/martina/Bigcat/Papers/2014-PLOSBiology-DiabeticLiver/analysis/network-analysis/drug-extension/drugbank.xml");
		try {
			System.out.println("[INFO:]\tRead DrugBank XML file");
			Set<DrugModel> drugs = new DrugBankParser().parse(drugBankXml);

			System.out.println("[INFO:]\t" + drugs.size() + " drugs are loaded.");
			
			System.out.println("[INFO:]\tSet up identifier mapping.");
			File bridgeFile = new File("/home/martina/Data/BridgeDb/Hs_Derby_20130701.bridge");
			Class.forName("org.bridgedb.rdb.IDMapperRdb");  
			IDMapper mapper = BridgeDb.connect("idmapper-pgdb:" + bridgeFile.getAbsolutePath());
			System.out.println("[INFO:]\tIdentifier mapping.database loaded");
			
			System.out.println("[INFO:]\tConvert DrugBank to XGMML graph.");
			Graph graph = convert(drugs, mapper);
			writeGraph(graph, new File("/home/martina/Bigcat/Papers/2014-PLOSBiology-DiabeticLiver/analysis/network-analysis/drug-extension/drugbank4.xgmml"));
			System.out.println("[INFO:]\tConversion is finished with " + graph.getNodes().size() + " nodes and " + graph.getEdges().size() + " edges.");
			
		} catch (JDOMException e) {
			System.out.println("[ERROR]\tCould not read drug bank file.");
		} catch (IOException e) {
			System.out.println("[ERROR]\tCould not find drug bank file.\t" + drugBankXml.getAbsolutePath());
		} catch (ClassNotFoundException e) {
			System.out.println("[ERROR]\tCould not set up identifier mapping.");
		} catch (IDMapperException e) {
			System.out.println("[ERROR]\tCould not set up identifier mapping.");
		}
		
	}
	
	private static Set<String> edges = new HashSet<String>();
	
	private static Graph convert(Set<DrugModel> drugs, IDMapper mapper) throws IDMapperException {
		Graph graph = new Graph();
		graph.setTitle("DrugBank_v4 (approved)");
		graph.setAttribute(CommonAttributes.DATABASE.getName(), "DrugBank_v4 (approved)");
		graph.setAttribute(CommonAttributes.TYPE.getName(), "drug-target interactions");
		
		Map<String, Node> nodes = new HashMap<String, Node>();
		
		for(DrugModel model : drugs) {
			if(model.getGroups().contains("approved")) {
				Node drug = graph.addNode(model.getDrugbankID());
				String identifiers = "[" + model.getDrugbankID();
				if(!model.getCasNumber().equals("")) identifiers = identifiers + "," + model.getCasNumber();
				identifiers = identifiers + "]";
				drug.appendAttribute("identifiers", identifiers);
				drug.appendAttribute("label", model.getName());
				drug.appendAttribute("name", model.getName());
				drug.appendAttribute("biologicalType", "drug");
				drug.appendAttribute("cas-number", model.getCasNumber());
				drug.appendAttribute("inchikey", model.getInChiKey());
				drug.appendAttribute("categories", model.getCategories().toString());
				drug.appendAttribute("groups", model.getGroups().toString());
				drug.appendAttribute("drugbank", model.getDrugbankID());
				
				if(!nodes.containsKey(model.getDrugbankID())) {
					nodes.put(model.getDrugbankID(), drug);
				} else {
					System.out.println("ERROR! Multiple nodes with same id!");
				}
				
				for(TargetModel target : model.getTargets()) {
					String uniprot = target.getUniprotId();
					String ensembl = "";
					if(!uniprot.equals("")) {
						Set<Xref> res = mapper.mapID(new Xref(uniprot, DataSource.getBySystemCode("S")), DataSource.getBySystemCode("En"));
						if(res.size() > 0) {
							ensembl = res.iterator().next().getId();
						}
					}
					if(!ensembl.equals("")) {
						Node gene;
						if(nodes.containsKey(ensembl)) {
							gene = nodes.get(ensembl);
						} else {
							gene = graph.addNode(ensembl);
							gene.appendAttribute("geneName", target.getGeneName());
							gene.appendAttribute("proteinName", target.getName());
							gene.appendAttribute("organism", target.getOrganism());
							gene.appendAttribute("uniProt", uniprot);
							gene.appendAttribute("ensembl", ensembl);
							gene.appendAttribute("biologicalType", "gene");
							
							String tIds = "[" + ensembl;
							Set<Xref> res = mapper.mapID(new Xref(ensembl, DataSource.getBySystemCode("En")), DataSource.getBySystemCode("L"));
							for(Xref x : res) {
								tIds = tIds + "," + x.getId();
							}
							Set<Xref> res2 = mapper.mapID(new Xref(ensembl, DataSource.getBySystemCode("En")), DataSource.getBySystemCode("S"));
							for(Xref x : res2) {
								tIds = tIds + "," + x.getId();
							}
							tIds = tIds + "]";
							gene.appendAttribute("identifiers", tIds);
						}
						nodes.put(ensembl, gene);
						
						addEdge(drug, gene, graph);
						
					}
					
					
				}
			}
		}
		
		return graph;
	}
	
	private static void addEdge(Node drug, Node gene, Graph graph) {
		String id = drug.getId() + "-" + gene.getId();
		if(!edges.contains(id)) {
			Edge e = graph.addEdge(id, drug, gene);
			e.setAttribute("datasource", "DrugBank_v4 (approved)");
			e.setAttribute("interactionType", "drug-target");
		}
	}
	
	private static void writeGraph(Graph graph, File output) throws IOException {
		PrintWriter po = new PrintWriter(output);
		XGMMLWriter.write(graph, po);
		po.close();
	}

}
