<template>
  <v-col>
    <v-toolbar dense class="toolbar">
      <v-toolbar-title>Users</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn v-if="isAdmin" text to="/main/admin/users/create" color="primary">Create User</v-btn>
      </v-toolbar-items>
    </v-toolbar>

    <v-card>
      <v-card-title>
        <v-spacer />
        <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="items"
        :loading="!items"
        :search="search"
        :items-per-page="15"
        :footer-props="{
          itemsPerPageOptions: [10, 15, 20, -1],
          showFirstLastPage: true,
          showCurrentPage: true,
        }"
        multi-sort
      >
        <template v-slot:item.is_active="{ item }">
          <v-icon v-if="item.is_active">mdi-check</v-icon>
        </template>
        <template v-slot:item.is_admin="{ item }">
          <v-icon v-if="item.is_admin">mdi-check</v-icon>
        </template>
        <template v-slot:item.action="{ item }">
          <v-menu bottom left>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on">
                <v-icon>mdi-dots-vertical</v-icon>
              </v-btn>
            </template>
            <v-list dense>
              <v-list-item
                :to="{
                  name: 'admin-users-edit',
                  params: {
                    id: item.id,
                  },
                }"
              >
                <v-list-item-icon>
                  <v-icon color="grey">mdi-pencil-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Edit</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list>
          </v-menu>
        </template>
      </v-data-table>
    </v-card>
  </v-col>
</template>

<script lang="ts">
import { userModule } from "@/modules/user";
import { Component, Vue } from "vue-property-decorator";
import { mainModule } from "@/modules/main";

@Component
export default class AdminUsers extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly userContext = userModule.context(this.$store);

  readonly headers = [
    {
      text: "Email",
      sortable: true,
      value: "email",
      align: "left",
    },
    {
      text: "Name",
      sortable: true,
      value: "name",
      align: "left",
    },
    {
      text: "Is Active",
      sortable: true,
      value: "is_active",
      align: "left",
    },
    {
      text: "Is Admin",
      sortable: true,
      value: "is_admin",
      align: "left",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
    },
  ];

  search = "";

  get isAdmin() {
    return this.mainContext.getters.isAdmin;
  }

  get items() {
    return this.userContext.getters.users;
  }

  async mounted() {
    await this.userContext.actions.getUsers();
  }
}
</script>

<style scoped>
.toolbar {
  margin-bottom: 10px;
}
</style>
