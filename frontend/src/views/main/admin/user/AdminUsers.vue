<template>
  <div>
    <v-toolbar dense>
      <v-toolbar-title>
        Manage Users
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn
        small
        color="primary"
        to="/main/admin/users/create"
      >Create User
      </v-btn>
    </v-toolbar>
    <v-data-table
      :headers="headers"
      :items="users"
    >
      <template v-slot:item.is_active="{ item }">
        <v-icon v-if="item.is_active">mdi-check</v-icon>
      </template>

      <template v-slot:item.is_superuser="{ item }">
        <v-icon v-if="item.is_superuser">mdi-check</v-icon>
      </template>

      <template v-slot:item.action="{ item }">
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon :to="{name: 'main-admin-users-edit', params: {id: item.id}}">
              <v-icon>mdi-pencil</v-icon>
            </v-btn>
          </template>
          <span>Edit</span>
        </v-tooltip>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
  import { userModule } from '@/modules/user';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class AdminUsers extends Vue {
    readonly userContext = userModule.context(this.$store);

    headers = [
      {
        text: 'Email',
        sortable: true,
        value: 'email',
        align: 'left',
      },
      {
        text: 'Full Name',
        sortable: true,
        value: 'full_name',
        align: 'left',
      },
      {
        text: 'Is Active',
        sortable: true,
        value: 'is_active',
        align: 'left',
      },
      {
        text: 'Is Superuser',
        sortable: true,
        value: 'is_superuser',
        align: 'left',
      },
      {
        text: 'Actions',
        value: 'action',
        sortable: false,
      },
    ];

    get users() {
      return this.userContext.getters.users;
    }

    async mounted() {
      await this.userContext.actions.getUsers();
    }
  }
</script>
