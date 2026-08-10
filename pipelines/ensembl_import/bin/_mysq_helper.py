"""
    MySQL helper functions for fetching and loading data from a MySQL database.
"""

from __future__ import annotations

import logging
from typing import Any, Sequence

import pymysql

logger = logging.getLogger(__name__)


def mysql_fetch_data(
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str = "",
    params: Sequence[Any] | None = None,
) -> list[tuple[Any, ...]]:
    """
    Run a simple SELECT query and return fetched rows.
    Returns an empty list on error.

    Args:
        query (str): SQL SELECT query to execute.
        database (str): Name of the database.
        host (str): Host address of the MySQL server.
        port (int): Port number of the MySQL server.
        user (str): Username to connect to the database.

    Returns:
        List of tuples representing query results.
    """
    conn = None
    cursor = None
    info: list[tuple[Any, ...]] = []
    try:
        conn = pymysql.connect(host=host, user=user, password=password, port=port, database=database.strip())
        cursor = conn.cursor()
        cursor.execute(query, params or ())
        info = list(cursor.fetchall())
    except pymysql.Error:
        logger.exception("MySQL error while executing query: %s", query)
    finally:
        if cursor is not None:
            try:
                cursor.close()
            except Exception:
                logger.exception("Error closing cursor")
        if conn is not None:
            try:
                conn.close()
            except Exception:
                logger.exception("Error closing connection")
    return info


def mysql_execute_query(
    query: str,
    database: str,
    host: str,
    port: int,
    user: str,
    password: str = "",
    params: Sequence[Any] | None = None,
) -> None:
    """
    Execute a SQL query without returning any results.
    """
    conn = None
    cursor = None
    try:
        conn = pymysql.connect(host=host, user=user, password=password, port=port, database=database.strip())
        cursor = conn.cursor()
        cursor.execute(query, params or ())
        conn.commit()
    except pymysql.Error:
        logger.exception("MySQL error while executing query: %s", query)
        raise
    finally:
        if cursor is not None:
            try:
                cursor.close()
            except Exception:
                logger.exception("Error closing cursor")
        if conn is not None:
            try:
                conn.close()
            except Exception:
                logger.exception("Error closing connection")
