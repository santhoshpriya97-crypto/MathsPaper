import re, time, networkx as nx

def closed_neighborhood_bit(adj, v):
    mask = 1 << v
    for u, conn in enumerate(adj[v]):
        if conn:
            mask |= 1 << u
    return mask

def dominates_bit(adj, D_mask, n):
    dom = 0
    for v in range(n):
        if (D_mask >> v) & 1:
            dom |= closed_neighborhood_bit(adj, v)
    return dom == (1 << n) - 1

def count_enclaves_bit(adj, D_mask, n):
    cnt = 0
    for v in range(n):
        if (D_mask >> v) & 1 and (closed_neighborhood_bit(adj, v) & ~D_mask) == 0:
            cnt += 1
    return cnt

def find_minimum_enclave_sets_one_per_enclave(adj):
    n = len(adj)
    if n <= 1:
        return (n, [0], [1]) if n else (0, [], [])
    closed = [closed_neighborhood_bit(adj, v) for v in range(n)]
    best_for = {}
    full = (1 << n) - 1

    def search(cur, enc):
        sz = cur.bit_count()
        if enc in best_for and sz >= best_for[enc][0]:
            return
        if dominates_bit(adj, cur, n) and count_enclaves_bit(adj, cur, n) == 1:
            best_for[enc] = (sz, cur)
            return

        dom = 0
        for v in range(n):
            if (cur >> v) & 1:
                dom |= closed[v]
        undom = full ^ dom
        if not undom:
            return
        cand = [v for v in range(n) if not ((cur >> v) & 1)]
        cand.sort(key=lambda x: -bin(closed[x] & undom).count('1'))
        for v in cand:
            new = cur | (1 << v)
            if count_enclaves_bit(adj, new, n) != 1:
                continue
            search(new, enc)

    for u in range(n):
        base = closed[u]
        if count_enclaves_bit(adj, base, n) == 1:
            search(base, u)

    if not best_for:
        return 0, [], []
    gmin = min(best_for[u][0] for u in best_for)
    verts = [u for u in best_for if best_for[u][0] == gmin]
    masks = [best_for[u][1] for u in verts]
    return gmin, sorted(verts), masks

def output_results(adj, graph_label="G"):
    start_time = time.perf_counter()
    gamma, verts, masks = find_minimum_enclave_sets_one_per_enclave(adj)
    end_time = time.perf_counter()
    elapsed = end_time - start_time
    print(f"Running time for {graph_label}: {elapsed:.6f} seconds")
    if gamma == 0:
        print(f"γ_ε({graph_label}) does not exist")
        return
    n = len(adj)
    def label(i):
        return f"u_{i+1}"
    print(f"γ_ε({graph_label}) = {gamma}")
    print("Enclave dominating vertices:", ", ".join(label(v) for v in verts))
    print(f"Total number of enclave dominating vertices: {len(verts)}")
    print("Minimum enclave dominating sets:")
    for m in masks:
        elems = sorted(label(v) for v in range(n) if (m >> v) & 1)
        print("{" + ", ".join(elems) + "}")

def graph_from_name(name):
    name = name.strip().upper()
    if name.startswith('P'):
        G = nx.path_graph(int(name[1:]))
    elif name.startswith('C'):
        G = nx.cycle_graph(int(name[1:]))
    elif name.startswith('K') and ',' not in name:
        G = nx.complete_graph(int(name[1:]))
    elif name.startswith('S'):
        n = int(name[1:])
        G = nx.star_graph(n - 1)
    elif name.startswith('W'):
        G = nx.wheel_graph(int(name[1:]))
    elif name.startswith('K') and ',' in name:
        a, b = map(int, name[1:].split(','))
        G = nx.complete_bipartite_graph(a, b)
    else:
        raise ValueError("Unknown graph name")
    n = G.number_of_nodes()
    adj = [[0] * n for _ in range(n)]
    for u, v in G.edges():
        adj[u][v] = adj[v][u] = 1
    return adj

def build_adjacency_matrix(inp):
    if isinstance(inp, list) and all(isinstance(r,list) for r in inp):
        return inp
    if isinstance(inp, tuple) and len(inp)==2:
        n, edges = inp
        mat = [[0]*n for _ in range(n)]
        for u,v in edges:
            mat[u][v]=mat[v][u]=1
        return mat
    if isinstance(inp, str):
        if re.match(r'^[PCKSW]\d+$|^K\d+,\d+$', inp):
            return graph_from_name(inp)
        parts = inp.strip().split()
        n = int(parts[0])
        mat = [[0]*n for _ in range(n)]
        for tok in parts[1:]:
            u,v = map(int, tok.split('-'))
            mat[u-1][v-1]=mat[v-1][u-1]=1
        return mat
    raise ValueError("Unsupported input")

def enclave_domination_from_input(inp):
    if isinstance(inp, str) and re.match(r'^[PCKSW]\d+$|^K\d+,\d+$', inp):
        label = inp.upper()
    elif isinstance(inp, str):
        label = "custom"
    elif isinstance(inp, tuple):
        label = f"graph_{inp[0]}"
    elif isinstance(inp, list):
        label = "adj_matrix"
    else:
        label = "G"
    adj = build_adjacency_matrix(inp)
    output_results(adj, label)

if __name__ == "__main__":
    enclave_domination_from_input("P15")
