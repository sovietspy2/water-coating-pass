echo "Creating .db file"
sqlite3 database.db ".import --csv ../research_history.csv runs"
echo "Loading .db file to temporary SQL db and starting sqlite3, you can find the data is the runs table!"
sqlite3 database.db