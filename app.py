# run the file on terminal on code location
# py -m streamlit run "path\to\file\app.py"
########################

import streamlit as st
from pymongo import MongoClient
from neo4j import GraphDatabase
from graphdatascience import GraphDataScience
from pyvis.network import Network

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# MongoDB Connection
client = MongoClient('mongodb://localhost:27017/')
db = client['Project-1st-Semester']
collection = db['protein_data']

# Checking whether is working or not
# print(collection.count_documents({}))

# Neo4j Connection
driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "rezareza"))

#################################################################################################################


# MongoDB Queries
def search_protein_mongodb(entry=None, interpro=None, ec_number=None):
    query = {}
    
    # Build query based on input parameters
    if entry:
        query['Entry'] = entry  # Exact match for Entry
    if interpro:
        query['InterPro'] = {'$regex': interpro, '$options': 'i'}  # Case-insensitive substring match for InterPro
    if ec_number:
        query['EC number'] = {'$regex': ec_number, '$options': 'i'}  # Case-insensitive substring match for EC number

    # Perform the query in MongoDB
    results = collection.find(query)
    return list(results)


# MongoDB Queries - Counts
def count_protein_mongodb(entry=None, interpro=None, ec_number=None):
    query = {}
    
    # Build query based on input parameters
    if entry:
        query['Entry'] = entry  # Exact match for Entry
    if interpro:
        query['InterPro'] = {'$regex': interpro, '$options': 'i'}  # Case-insensitive substring match for InterPro
    if ec_number:
        query['EC number'] = {'$regex': ec_number, '$options': 'i'}  # Case-insensitive substring match for EC number
    
    # Perform the query in MongoDB
    results = collection.count_documents(query)
    return results 

# Neo4j Function to Search Proteins
def search_protein_neo4j(entry=None, interpro=None, ec_number=None):
    query = """
    MATCH (p:Protein1)
    WHERE 
        ($entry IS NULL OR COALESCE(p.id, "") CONTAINS $entry) AND
        ($interpro IS NULL OR COALESCE(p.interpro, "") CONTAINS $interpro) AND
        ($ec_number IS NULL OR COALESCE(p.ec_number, "") CONTAINS $ec_number)
    RETURN 
        COALESCE(p.id, "N/A") AS Entry,
        COALESCE(p.interpro, "N/A") AS InterPro,
        COALESCE(p.ec_number, "N/A") AS ECNumber
    """
    with driver.session() as session:
        result = session.run(query, entry=entry, interpro=interpro, ec_number=ec_number)
        return [{"Entry": record["Entry"], "InterPro": record["InterPro"], "ECNumber": record["ECNumber"]} for record in result]



# Neo4j Function to Retrieve Neighbors
def get_protein_neighbors(entry):
    query = """
    MATCH (p:Protein1 {id: $entry})-[r:DOMAIN_SIMILARITY]-(neighbor:Protein1)
    RETURN neighbor.id AS NeighborID, r.weight AS Weight
    """
    with driver.session() as session:
        result = session.run(query, entry=entry)
        neighbors = [{"Neighbor": record["NeighborID"], "Weight": record["Weight"]} for record in result]
        return neighbors



# Function to build and render a protein graph
def build_protein_graph(entry):
    query = """
    MATCH (p:Protein1 {id: $entry})-[r:DOMAIN_SIMILARITY]-(neighbor:Protein1)
    RETURN p.id AS Source, neighbor.id AS Target, r.weight AS Weight
    """
    with driver.session() as session:
        result = session.run(query, entry=entry)

        # Convert query result into a list of dictionaries
        graph_data = [{"Source": record["Source"], "Target": record["Target"], "Weight": record["Weight"]}
                      for record in result]

    # Initialize the Pyvis Network graph
    net = Network(height="600px", width="100%", notebook=False, directed=False)

    # Add nodes and edges based on Neo4j query result
    nodes = set()
    core_protein = entry  # The queried core protein

    for edge in graph_data:
        source = edge["Source"]
        target = edge["Target"]
        weight = edge["Weight"]

        # Add core protein in green with black text inside
        if source == core_protein and source not in nodes:
            net.add_node(source, label=source, title=f"Protein: {source}", color="green", font={"size": 15, "color": "black"})
            nodes.add(source)

        # Add neighbors in red with black text inside
        if target not in nodes:
            net.add_node(target, label=target, title=f"Protein: {target}", color="red", font={"size": 15, "color": "black"})
            nodes.add(target)

        # Add edges with weights displayed as labels
        net.add_edge(source, target, value=weight, title=f"Weight: {weight:.4f}", color="blue", width=weight * 1.2)
        # net.add_edge(source, target, value=weight, labels=f"Weight: {weight}", color="blue", width=weight * 2)

    # Save graph to an HTML file
    net.save_graph("protein_graph.html")
    return "protein_graph.html"



# Function to build a limited protein graph
def build_protein_graph_number(limit):
    query = """
    MATCH p=(n)-[r:DOMAIN_SIMILARITY]-(m)
    RETURN p, r.weight AS Weight
    LIMIT $limit
    """
    with driver.session() as session:
        result = session.run(query, limit=limit)

        # Extract nodes and relationships
        graph_data = []
        for record in result:
            for segment in record["p"].relationships:  # Each relationship in the path
                graph_data.append({
                    "Source": segment.start_node["id"],
                    "Target": segment.end_node["id"],
                    "Weight": record["Weight"]
                })

    # Initialize the Pyvis Network graph
    net = Network(height="600px", width="100%", notebook=False, directed=False)

    # Add nodes and edges
    nodes = set()
    for edge in graph_data:
        source = edge["Source"]
        target = edge["Target"]
        weight = edge["Weight"]

        # Add nodes (all in red with black text)
        if source not in nodes:
            net.add_node(source, label=source, title=f"Protein: {source}", color="red", font={"size": 15, "color": "black"})
            nodes.add(source)
        if target not in nodes:
            net.add_node(target, label=target, title=f"Protein: {target}", color="red", font={"size": 15, "color": "black"})
            nodes.add(target)

        # Add edges with weight labels
        net.add_edge(source, target, value=weight, title=f"Weight: {weight:.4f}", color="blue", width=weight * 1.2)

    # Save graph to an HTML file
    net.save_graph("limited_protein_graph.html")
    return "limited_protein_graph.html"




# Neo4j Function to Compute EC Statistics
def compute_protein_ec_statistics():
    query = """
    MATCH (p:Protein1)
    WITH p,
         CASE 
             WHEN p.ec_number IS NULL THEN 'No EC number'
             WHEN SIZE(split(p.ec_number, ';')) = 1 THEN 'One EC number'
             ELSE 'Multiple EC numbers'
         END AS ec_category
    RETURN ec_category, COUNT(*) AS count
    """
    with driver.session() as session:
        result = session.run(query)
        stats = {record["ec_category"]: record["count"] for record in result}
        return stats


# Function to fetch sequence lengths from MongoDB
def get_sequence_lengths():
    # Fetch all sequences from MongoDB
    sequences = collection.find({}, {"Sequence": 1, "_id": 0})
    
    # Calculate sequence lengths
    sequence_lengths = [len(doc["Sequence"]) for doc in sequences if "Sequence" in doc and doc["Sequence"]]
    
    return sequence_lengths



########################################################################################################################################

# Streamlit GUI
st.title("Final Project - Protein Database")

# MongoDB
st.sidebar.header("MongoDB")

# MongoDB Search
st.sidebar.subheader("Search MongoDB")
entry = st.sidebar.text_input("Entry", key="search_entry")  # Unique key
interpro = st.sidebar.text_input("InterPro", key="search_interpro")  # Unique key
ec_number = st.sidebar.text_input("EC Number", key="search_ec_number")  # Unique key

if st.sidebar.button("Search MongoDB"):
    results = search_protein_mongodb(entry=entry, interpro=interpro, ec_number=ec_number)
    if results:
        st.subheader("MongoDB Search Results")
        for result in results:
            st.json(result)  # Display results in JSON format for clarity
    else:
        st.write("No results found.")


# MongoDB Count
st.sidebar.subheader("Count MongoDB")
entry = st.sidebar.text_input("Entry", key="count_entry")  # Unique key
interpro = st.sidebar.text_input("InterPro", key="count_interpro")  # Unique key
ec_number = st.sidebar.text_input("EC Number", key="count_ec_number")  # Unique key

if st.sidebar.button("Count MongoDB"):
    result = count_protein_mongodb(entry=entry, interpro=interpro, ec_number=ec_number)
    if result:
        st.subheader("MongoDB Count Results")
        st.write(f"Number of matching documents: {result}")
    else:
        st.write("No results found.")



# Neo4j
st.sidebar.header("Neo4j")

# Neo4j Search
st.sidebar.subheader("Search Neo4j")
neo4j_entry = st.sidebar.text_input("Entry", key="neo4j_search_entry")
neo4j_interpro = st.sidebar.text_input("InterPro", key="neo4j_search_interpro")
neo4j_ec_number = st.sidebar.text_input("EC Number", key="neo4j_search_ec_number")

if st.sidebar.button("Search Neo4j"):
    results = search_protein_neo4j(entry=neo4j_entry, interpro=neo4j_interpro, ec_number=neo4j_ec_number)
    if results:
        st.subheader("Neo4j Search Results")
        for record in results:
            st.json(record)  # Display results in JSON format with meaningful keys
    else:
        st.write("No results found.")



# Neo4j Neighbor Retrieval
st.sidebar.subheader("Neighbors of a Protein")
neighbor_entry = st.sidebar.text_input("Protein Entry for Neighbors")

if st.sidebar.button("Show Neighbors"):
    if neighbor_entry:
        st.subheader(f"Neighbors for Protein: {neighbor_entry}")
        neighbors = get_protein_neighbors(neighbor_entry)
        if neighbors:
            for neighbor in neighbors:
                st.write(f"Neighbor:\t\t {neighbor['Neighbor']}, Weight:\t\t {neighbor['Weight']}")
        else:
            st.write("No neighbors found.")
    else:
        st.write("Please enter a valid protein entry.")


# Neo4j Graph Visualization
st.sidebar.subheader("Visualize Protein Graph")
graph_entry = st.sidebar.text_input("Protein Entry for Graph Visualization")

if st.sidebar.button("Show Protein Graph"):
    if graph_entry:
        st.subheader(f"Protein Graph for Entry: {graph_entry}")
        graph_file = build_protein_graph(graph_entry)
        st.components.v1.html(open(graph_file, "r", encoding="utf-8").read(), height=550)
    else:
        st.write("Please enter a valid protein entry.")



# Neo4j Graph Visualization Limit
st.sidebar.subheader("Visualize Limited Protein Graph")
graph_limit = st.sidebar.number_input("Number of Connections (Limit)", min_value=1, step=1)

if st.sidebar.button("Show Limited Protein Graph"):
    if graph_limit:
        st.subheader(f"Protein Graph with {graph_limit} Connections")
        graph_file = build_protein_graph_number(graph_limit)
        st.components.v1.html(open(graph_file, "r", encoding="utf-8").read(), height=550)
    else:
        st.write("Please enter a valid limit.")



# Statistics
st.sidebar.header("Statistics")

# Statistics 1: EC Number Distribution
if st.sidebar.button("EC Number Pie Chart - Sample Dataset"):
    st.subheader("Protein EC Number Distribution")
    
    # Fetch statistics from Neo4j
    ec_stats = compute_protein_ec_statistics()
    
    # Prepare labels and values for pie chart
    labels = list(ec_stats.keys())
    sizes = list(ec_stats.values())
    colors = ["#FF9999", "#66B3FF", "#99FF99"]  # Different colors for visualization
    
    # Create a pie chart
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=[f"{label} ({size})" for label, size in zip(labels, sizes)], autopct='%1.1f%%', startangle=140, colors=colors)

    ax.axis("equal")  # Equal aspect ratio ensures the pie is drawn as a circle
    
    # Display the pie chart
    st.pyplot(fig)


# Statistics 2: Sequence Length Distribution
if st.sidebar.button("Sequence Length Distribution - Main Dataset"):
    st.subheader("Sequence Length Distribution of Proteins")
    
    # Get sequence lengths from MongoDB
    sequence_lengths = get_sequence_lengths()
    
    if sequence_lengths:
        # Filter out zero or invalid lengths
        filtered_lengths = [x for x in sequence_lengths if x > 0]

        # Calculate bins dynamically using numpy
        min_length = min(filtered_lengths)
        max_length = max(filtered_lengths)
        bin_width = 200  # Set a reasonable bin width
        bins = np.arange(min_length, max_length + bin_width, bin_width)

        # Plot Barplot (Histogram) and KDE Curve
        fig, ax = plt.subplots(figsize=(10, 6))

        # Histogram with proper bins
        ax.hist(filtered_lengths, bins=bins, edgecolor='black', alpha=0.7, color="green", density=False, label="Bar Plot")

        # Set X-axis range
        ax.set_xlim(0, 4000)
        
        # Add labels, title, and legend
        ax.set_title("Protein Sequence Length Distribution")
        ax.set_xlabel("Sequence Length")
        ax.set_ylabel("Frequency")
        ax.legend()

        # Display the plot in Streamlit
        st.pyplot(fig)
    else:
        st.write("No sequence data found in the database.")