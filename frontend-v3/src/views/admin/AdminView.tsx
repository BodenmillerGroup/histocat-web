import { Route, Switch } from "react-router-dom";
import React from "react";
import { UsersView } from "./users/UsersView";
import { useProfileStore } from "../../modules/profile";
import { AccessDeniedView } from "../../components/AccessDeniedView";
import { ModelsView } from "./models/ModelsView";

export function AdminView() {
  const hasAdminAccess = useProfileStore((state) => state.hasAdminAccess());
  if (!hasAdminAccess) {
    return <AccessDeniedView title="Access Denied" content="You should have super-admin access rights." />;
  }
  return (
    <Switch>
      <Route path="/admin/users" component={UsersView} />
      <Route path="/admin/models" component={ModelsView} />
    </Switch>
  );
}
