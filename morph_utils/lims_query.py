
import psycopg2

CONNECTION_STRING = 'host=limsdb2 dbname=lims2 user=limsreader password=limsro'

def query(sql, args):
    """
    Query LIMS

    :param sql: SQL command
    :param args: params for SQL command
    :return: query result
    """
    conn = psycopg2.connect(CONNECTION_STRING)
    cur = conn.cursor()

    cur.execute(sql, args)
    results = cur.fetchall()

    cur.close()
    conn.close()

    return results


def get_name_by_id(specimen_id):
    """
    Get cell name from id

    :param specimen_id: cell specimen id
    :return: cell specimen name 
    """

    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.id = '{}'
    """.format(specimen_id)
    
    x = query(sql, ())[0][0]
    return x

def get_id_by_name(specimen_name):
    """
    Get cell id from name

    :param specimen_name: cell specimen name
    :return: cell specimen id 
    """

    sql = """
    SELECT sp.name as sp, sp.id
    FROM specimens sp
    WHERE sp.name = '{}'
    """.format(specimen_name)
    
    x = query(sql, ())[0][1]
    return x

#convert px -> um 
def query_for_z_resolution(specimen_id):
    """
    Get resolution of z slice axis for a cell. 

    :param specimen_id: cell specimen id
    :return: z axis resolution
    """
     
    sql = """
    select ss.id, ss.name, shs.thickness from specimens ss
    join specimen_blocks sb on ss.id = sb.specimen_id
    join blocks bs on bs.id = sb.block_id
    join thicknesses shs on shs.id = bs.thickness_id 
    where ss.id = {}
    """.format(specimen_id)
    
    res = query(sql,())
    try:
        return res[0][-1]
    except:
        return None
