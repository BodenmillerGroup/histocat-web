import { Redirect, Route, Switch } from "react-router-dom";
import React from "react";
import { useGroupsStore } from "modules/groups";
import { AccessDeniedView } from "components/AccessDeniedView";
import { GroupsListView } from "./GroupsListView";

export default function GroupsView() {
  // const isGroupAdmin = useGroupsStore((state) => state.isGroupAdmin());
  // if (!isGroupAdmin) {
  //   return <AccessDeniedView title="Access Denied" content="You should have super-admin access rights." />;
  // }
  return (
    <Switch>
      <Route path="/" component={GroupsListView} />
    </Switch>
  );
}
