
import psycopg2

CONNECTION_STRING = 'host=limsdb2 dbname=lims2 user=limsreader password=limsro'

def query(sql, args):
    conn = psycopg2.connect(CONNECTION_STRING)
    cur = conn.cursor()

    cur.execute(sql, args)
    results = cur.fetchall()

    cur.close()
    conn.close()

    return results


def get_name_by_id(sp_id):
    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.id = '{}'
    """.format(sp_id)
    
    x = query(sql, ())[0][0]
    return x

def get_id_by_name(sp_name):
    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.name = '{}'
    """.format(sp_name)
    
    x = query(sql, ())[0][1]
    return x
