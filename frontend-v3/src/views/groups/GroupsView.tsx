import { Route, Switch, useRouteMatch } from "react-router-dom";
import React, { useEffect, useState } from "react";
import { useGroupsStore } from "modules/groups";
import { GroupsListView } from "./GroupsListView";
import shallow from "zustand/shallow";
import { LoadingView } from "components/LoadingView";
import GroupView from "./GroupView";
import { WebSocketManager } from "utils/WebSocketManager";
import { useAuthStore } from "../../modules/auth";

export default function GroupsView() {
  const { path } = useRouteMatch();
  const token = useAuthStore(state => state.token);
  const { getGroups, getTags } = useGroupsStore(
    (state) => ({
      getGroups: state.getGroups,
      getTags: state.getTags,
    }),
    shallow
  );
  const [loading, setLoading] = useState<boolean>(true);

  // Mounted
  useEffect(() => {
    WebSocketManager.init(token!);
    Promise.all([getGroups(), getTags()]).then(() => setLoading(false));
  }, [getGroups, getTags, token]);

  // Unmounted
  useEffect(() => {
    return () => {
      WebSocketManager.close();
    };
  }, []);

  if (loading) {
    return <LoadingView />;
  }

  return (
    <Switch>
      <Route exact={true} path={path} component={GroupsListView} />
      <Route path={`${path}:groupId`} component={GroupView} />
    </Switch>
  );
}
